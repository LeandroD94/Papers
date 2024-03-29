# operating system: Ubuntu 22.04 (Linux)
# programs required: conda, qiime2-2022.8, bowtie-2.2.5, picrust-2.5, lefse-1.1.2, R-4.2

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

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 275 --p-trunc-len-r 168 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 7

qiime tools export --input-path denoising-stats.qza --output-path ./      ### to export and create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > Initial_and_processed_reads.tsv
cat Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9


chmod +rwx stats.tsv
rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv


################# REMOVING EUCARYOTIC ASV (CONTAMINANTS) #########################

qiime tools export --input-path rep-seqs.qza --output-path ./ # it will export "dna-sequences.fasta"

conda deactivate

conda activate bowtie2

bowtie2 -x ~/Scrivania/Tools/Database/INDEXED_GRCh38_analysis_set_BOWTIE2/GRCh38_noalt_as -U dna-sequences.fasta --al eucar_contaminant_ASV.fasta -S bowtie2_SAM_output.txt -f -p 4 2> percentuals_aligned_bowtie2.txt

conda deactivate


mkdir eucar_contaminants

mv rep-seqs.qza eucar_contaminants/rep-seqs_with_euc_contam.qza
#cp rep-seqs.qza eucar_contaminants/rep-seqs_with_euc_contam.qza

grep -Fwv -f eucar_contaminant_ASV.fasta dna-sequences.fasta > decontaminated_seq.fasta

echo -e " Number of human contaminant ASV sequences:" $(grep ">" eucar_contaminant_ASV.fasta | wc -l ) "among" $(grep ">" decontaminated_seq.fasta | wc -l) "procaryotic ASVs" > eucar_contaminants/Number_of_eucar_Contaminant.txt

chmod +rwx dna-sequences.fasta
rm dna-sequences.fasta


mv eucar_contaminant_ASV.fasta eucar_contaminants/eucar_contaminant_ASV.fasta
mv bowtie2_SAM_output.txt eucar_contaminants/bowtie2_SAM_output.txt
mv percentuals_aligned_bowtie2.txt eucar_contaminants/percentuals_aligned_bowtie2.txt


conda activate qiime2-2022.8

qiime tools import --input-path decontaminated_seq.fasta --output-path clean_rep-seqs.qza --type FeatureData[Sequence]



######################### TAXONOMIC CLASSIFICATION ##################################       requires ~20GB RAM

### BAYESIAN CLASSIFICATION (using Scikit-learn retrained on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier ~/Scrivania/Tools/Database/SILVA_138_V3V4_BayesClassifier.qza --i-reads clean_rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization confidence_taxonomy.qzv



### BAYESIAN CLASSIFICATION ALSO ON THE HOST CONTAMINATED OBJECT
qiime feature-classifier classify-sklearn --i-classifier ~/Scrivania/Tools/Database/SILVA_138_V3V4_BayesClassifier.qza --i-reads eucar_contaminants/rep-seqs_with_euc_contam.qza --o-classification eucar_contaminants/contam_taxonomy.qza



################# COMPUTING THE PHYLOGENETIC TREE DE NOVO #############################

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



########################## PICRUST2-LEFSE ##############################

#### PREPARING FILES FOR PICRUST 

mkdir PICRUST2_LEFSE

qiime tools export --input-path clean_rep-seqs.qza --output-path PICRUST2_LEFSE # it will export "dna-sequences.fasta"
qiime tools export --input-path table.qza --output-path PICRUST2_LEFSE # it will export "feature-table.biom"

conda deactivate

cd PICRUST2_LEFSE


##### PICRUST2

conda activate picrust2

picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2 -p 4 --verbose -e 0.0 -t sepp 

# the last two argument after verbose are for SEPP algorithm ("e" needed to avoid and error, see https://github.com/gavinmdouglas/q2-picrust2/issues/13)

add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m KO -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz

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
if(identical(Metadata$FASTQ_ID,colnames(a))){a<-rbind.data.frame(Metadata$Condition,a)} else {cat ("\nsomething gone wrong...\n\n")}
a<-cbind.data.frame(c("Condition",Descriptions$pathway),a)
colnames(a)[colnames(a)=='c("Condition", Descriptions$pathway)']<-"Sample"
head(a,n=3)

write.table(a,file="Output_for_LefSe.tsv",row.names = F,quote = F, sep="\t")

rm(a, Descriptions, Metadata)

quit(save="no")



######### LEFSE

conda activate lefse

lefse_format_input.py 'Output_for_LefSe.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000 # groups is second row, normalization to 1 milion (see Galaxy recommendation)

lefse_run.py formatted_table_for_lefse.in -y 1 Result_LEFSE.res # y is for multiclass analysis

rm formatted_table_for_lefse.in

conda deactivate

cd ..
