######################## PROCESSING RAW SALIVA FASTQ ########################

conda activate qiime2-2021.4

mkdir QIIME_SALIVA

cd QIIME_SALIVA

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_SALIVA_dir --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

rm Log_di_cutadapt.txt


##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"

qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 10 --p-trim-left-r 0 --p-trunc-len-f 268 --p-trunc-len-r 194 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 8

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza
 
### BAYESIAN CLASSIFICATION (using retrained RDP on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier ~/Desktop/SILVA_138_V3V4_BayesClassifier.qza --i-reads rep-seqs.qza --o-classification Bayes_saliva_taxonomy.qza

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


############################# PICRUST2-LEFSE (SALIVA) ##############################

cd PICRUST2_LEFSE

conda activate picrust2

picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2 -p 8

add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz

# --> output labels modified in R (see analysis script) --> LefSe (Galaxy)

conda deactivate


##################### RESETTING ORIGINAL WORKING DIRECTORY ####################

cd ../..

######################## PROCESSING RAW STOOL FASTQ #########################

conda activate qiime2-2021.4

mkdir QIIME2_STOOL

cd QIIME2_STOOL

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_STOOL_dir --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

rm Log_di_cutadapt.txt


##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"

qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 261 --p-trunc-len-r 185 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 8

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza
 
### BAYESIAN CLASSIFICATION (using retrained RDP on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier ~/Desktop/SILVA_138_V3V4_BayesClassifier.qza --i-reads rep-seqs.qza --o-classification Bayes_feci_taxonomy.qza

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


############################# PICRUST2-LEFSE (STOOL) ##############################

cd PICRUST2_LEFSE

conda activate picrust2

picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2 -p 8

add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz

# --> output labels modified in R (see analysis script) --> LefSe (Galaxy)

conda deactivate
