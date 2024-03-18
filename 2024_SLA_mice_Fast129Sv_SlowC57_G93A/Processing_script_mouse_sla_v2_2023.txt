# Operative Sistem: Ubuntu 22.04 (Linux)
# required: conda, qiime2-2022.8, bowtie-2.2.5, picrust-2.5, lefse-1.1.2, R-4.2


########################## PREPARING THE ENVINRONMENT ############################

mkdir QIIME
chmod +rwx QIIME

##### FASTQ CHECKSUM

md5sum raw_FASTQ_dir/* | sed 's+raw_FASTQ_dir/++g' | sed 's/  / \t/' > QIIME/checksum_FASTQ.tsv

dir raw_FASTQ_dir/ | sed '/R2_001/d' >R1.txt
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt


############################ PROCESSING RAW FILES ##########################

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


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 274 --p-trunc-len-r 169 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 8

qiime tools export --input-path denoising-stats.qza --output-path ./      ### to export and create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > Initial_and_processed_reads.tsv
cat Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9

chmod +rwx stats.tsv


################# REMOVING EUCARYOTIC ASV (CONTAMINANTS) #########################

qiime tools export --input-path rep-seqs.qza --output-path ./ # it will export "dna-sequences.fasta"

conda deactivate

conda activate bowtie2

bowtie2 -x ~/Scrivania/Tools/Database/INDEXED_mouse_GRCm39_BOWTIE2/GRCm39 -U dna-sequences.fasta --al eucar_contaminant_ASV.fasta -S bowtie2_SAM_output.txt -f -p 4 2> percentuals_aligned_bowtie2.txt

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



######################### TAXONOMIC CLASSIFICATION ##################################


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




#################### PICRUST2 (FUNCTIONALITY PREDICTION) ########################


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

## TO CONVERT KEGG IN KO
pathway_pipeline.py -i picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -o picrust2/KEGG_pathways --map /home/leandro/miniconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv --no_regroup

add_descriptions.py -i picrust2/KEGG_pathways/path_abun_unstrat.tsv.gz -o picrust2/KEGG_pathways/path_abun_unstrat_descriptions.tsv.gz  --custom_map_table /home/leandro/miniconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz


conda deactivate



##### SUBSETTING FEATURES FOR LEFSE

R # open R

a <- read.delim("picrust2/KEGG_pathways/path_abun_unstrat_descriptions.tsv.gz")

if (!require('stringr')) {install.packages('stringr')}

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]

### Names update
sample<-substring(colnames(a), first = 10, last= 15) # strtrim to cut until X
sample<-gsub(".A","",sample)
sample<-gsub(".P","",sample)
sample<-gsub(".","",sample, fixed=T)
colnames(a)<-sample


### Importing infos
Metadata <- as.data.frame(read.table("../../Metadata.tsv", sep="\t", header=T))
#head(Metadata)
Metadata$FASTQ_ID<-as.character(Metadata$FASTQ_ID)
rownames(Metadata)<-Metadata$FASTQ_ID
Metadata<-Metadata[colnames(a),]
if(identical(Metadata$FASTQ,colnames(a))){a<-rbind.data.frame(Metadata$Strain,Metadata$Genotype,a)} else {stop ("\nsomething gone wrong...\n\n")}
a<-cbind.data.frame(c("Strain","Genotype",Descriptions$pathway),a)
colnames(a)[colnames(a)=='c("Strain", "Genotype", Descriptions$pathway)']<-"Sample"
#head(a,n=3)


### Subsetting
PS<-row.names(Metadata[ Metadata$Time_point=="PS", ])
G93A<-row.names(Metadata[ Metadata$Strain=="G93A", ])   # vs Ntg
C57Ola <- row.names(Metadata[ Metadata$Genotype=="C57Ola", ])  # vs 129Sv

# PS G93A --> 129sv vs C57Ola
targets<- colnames(a)[ colnames(a) %in% PS & colnames(a) %in% G93A ]
table_target <- a[ a[,1]!="Strain" , ]
write.table( table_target[ , c("Sample",targets)], file="Output_for_LefSe_inTimePS_inStrainG93A_C57Ola_vs_129Sv.tsv",row.names = F,quote = F, sep="\t")

# PS Ntg --> 129sv vs C57Ola
targets<- colnames(a)[ colnames(a) %in% PS & !(colnames(a) %in% G93A) ]
table_target <- a[ a[,1]!="Strain" , ]
write.table( table_target[ , c("Sample",targets)], file="Output_for_LefSe_inTimePS_inStrainNtg_C57Ola_vs_129Sv.tsv",row.names = F,quote = F, sep="\t")

# PS C57Ola --> Ntg vs G93A
targets<- colnames(a)[ colnames(a) %in% PS & colnames(a) %in% C57Ola ]
table_target <- a[ a[,1]!="Genotype" , ]
write.table( table_target[ , c("Sample",targets)], file="Output_for_LefSe_inTimePS_inGenotypeC57Ola_Ntg_vs_G93A.tsv",row.names = F,quote = F, sep="\t")

# PS 129Sv --> Ntg vs G93A
targets<- colnames(a)[ colnames(a) %in% PS & !(colnames(a) %in% C57Ola) ]
table_target <- a[ a[,1]!="Genotype" , ]
write.table( table_target[ , c("Sample",targets)], file="Output_for_LefSe_inTimePS_inGenotype129Sv_Ntg_vs_G93A.tsv",row.names = F,quote = F, sep="\t")


q("no")



##### STATISTICAL ANALYSES WITH LEFSE

conda activate lefse

for i in Output_for_LefSe_inTimePS_* 
do
sample_name=$( echo ${i/*_inTime/} ) # NB: still with ".tsv" in it

echo $i 
echo -e "\n"
lefse_format_input.py  $i  formatted_table_for_lefse.in -c 2 -u 1 -o 1000000 # groups is second row, normalization to 1 milion (see Galaxy recommendation)

lefse_run.py formatted_table_for_lefse.in -y 1 Result_LEFSE_${sample_name/.tsv/}.res # y is for multiclass analysis

rm formatted_table_for_lefse.in

done



cd ..   # PICRUST2 step concluded, returning to the QIIME folder



############# EXPORTING ORIGINAL NUMBER OF READS ############################

qiime demux summarize --i-data demux-paired-end.qza --o-visualization original_reads.qzv 
qiime tools export --input-path original_reads.qzv  --output-path Raw_Reads_Info
mv Raw_Reads_Info/per-sample-fastq-counts.tsv Original_number_of_reads_for_sample.tsv
chmod +rw Original_number_of_reads_for_sample.tsv
chmod +rw Raw_Reads_Info -R
rm Raw_Reads_Info -R


############# CLEANING THE FOLDER FROM QIIME TEMP FILES #####################

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv original_reads.qzv 


