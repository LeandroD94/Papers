# Operative System: Ubuntu 22.04 (Linux)
# required: conda, qiime2-2022.8, bowtie-2.2.5, picrust-2.5, lefse-1.1.2, R-4.2

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

conda activate qiime2-2022.8

cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_dir --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza


##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --p-minimum-length 100 --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

cat Log_reads_with_primers.txt


##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 200000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 276 --p-trunc-len-r 196 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 7

qiime tools export --input-path denoising-stats.qza --output-path ./      ### to export and create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > DADA2_Initial_and_processed_reads.tsv
cat DADA2_Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9


chmod +rwx stats.tsv
rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv


######### CHECK OF HUMAN CONTAMINANTS WITH BOWTIE2 ...

# NO HUMAN CONTAMINANTS !

mv rep-seqs.qza clean_rep-seqs.qza




######################### TAXONOMIC CLASSIFICATION ##################################

### BAYESIAN CLASSIFICATION (using Scikit-learn retrained on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier ~/Scrivania/Tools/Database/SILVA_138_V3V4_BayesClassifier.qza --i-reads clean_rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization confidence_taxonomy.qzv




################# COMPUTING THE PHYLOGENETIC TREE DE NOVO #############################


#### DE NOVO TREE ALIGNMENT

qiime alignment mafft --i-sequences clean_rep-seqs.qza --o-alignment aligned-seqs.qza
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



########################## PICRUST2-LEFSE ##############################

#### PREPARING FILES FOR PICRUST 

mkdir PICRUST2

qiime tools export --input-path clean_rep-seqs.qza --output-path PICRUST2 # it will export "dna-sequences.fasta"
qiime tools export --input-path table.qza --output-path PICRUST2 # it will export "feature-table.biom"

conda deactivate

cd PICRUST2


##### PICRUST2

conda activate picrust2

picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2 -p 4 --verbose -e 0.0 -t sepp 
# the last two argument after verbose are for SEPP algorithm ("e" needed to avoid and error, see https://github.com/gavinmdouglas/q2-picrust2/issues/13)

# add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m metacyc  -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz


## TO CONVERT KEGG IN KO
pathway_pipeline.py -i picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -o picrust2/KEGG_pathways --map /home/leandro/miniconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv --no_regroup

add_descriptions.py -i picrust2/KEGG_pathways/path_abun_unstrat.tsv.gz -o picrust2/KEGG_pathways/path_abun_unstrat_descriptions.tsv.gz  --custom_map_table /home/leandro/miniconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz


conda deactivate



