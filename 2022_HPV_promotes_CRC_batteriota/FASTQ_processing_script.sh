# requires qiime2 running in conda environment through bash (Ubuntu 20.08)


mkdir QIIME

##### FASTQ CHECKSUM

md5sum raw_FASTQ_dir/* | sed 's+raw_FASTQ_dir/++g' | sed 's/  / \t/' > QIIME/checksum_FASTQ.tsv

dir raw_FASTQ_dir/ | sed '/R2_001/d' >R1.txt
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt

##### PROCESSING

conda activate qiime2-2021.4

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path raw_FASTQ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 8 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000

qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 261 --p-trunc-len-r 198 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 8

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv


qiime feature-classifier classify-sklearn --i-classifier SILVA_138_V3V4_BayesClassifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization confidence_taxonomy.qzv


