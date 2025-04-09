#!/user/bin/env bash

# Operative System: Ubuntu (Linux)
# required: conda, fastqc, multiqc, bbmap (bbduk.sh) 



####### SETTING VARIABLES


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows loads the environment variable "into the script"

### BBDUK ADAPTERS FILE
ADAPTERS_FILE=~/Desktop/programs/bbmap_v39/resources/adapters.fa




### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### THREADS TO USE
threads=62
# nproc # to check the PC availability





####### PREPARING THE DIRECTORIES #######

mkdir processed_FASTQ 
mkdir FASTQ_check
# mkdir FASTQ_check/fastqc_raw_reads
mkdir FASTQ_check/BBDuk
mkdir FASTQ_check/BBDuk/Discarded




####### INITIAL QUALITY CONTROL #######

#for i in ./raw_FASTQ_dir/*gz ; do fastqc $i -o FASTQ_check/fastqc_raw_reads/ --threads $threads ;  echo -e "\n"  ; done

#multiqc FASTQ_check/fastqc_raw_reads/ -o FASTQ_check/fastqc_raw_reads

#open FASTQ_check/fastqc_raw_reads/multiqc_report.html
#cat FASTQ_check/fastqc_raw_reads/multiqc_data/multiqc_general_stats.txt | cut -f 1,3,5,7,9 | sed 's/FastQC_mqc-generalstats-fastqc-//g' > #FASTQ_check/fastqc_raw_reads/General_info_Sequences_recap.txt

#cat FASTQ_check/fastqc_raw_reads/General_info_Sequences_recap.txt

#rm FASTQ_check/fastqc_raw_reads/*fastqc.zip



# fastqc does not work by remote (?) but I still need this info
for i in ./DNAseq_raw_fastq/*gz
do 
echo "$i : $(zcat $i | grep "^@" -c)" >> FASTQ_check/original_number_of_seqs.txt
done
sed -e "s#./DNAseq_raw_fastq/##g" -e "s/.gz//g" -i FASTQ_check/original_number_of_seqs.txt




####### REMOVING THE BAD SEQUENCES WITH BBDUK #######

FOLDER_INPUT="DNAseq_raw_fastq"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep ".fastq" | extract_name )


for i in $SAMPLE_LIST    ### beginning of the loop
do

echo -e "\n\n\n\n *** Preprocessing the sample $i with bbduk... *** \n"

### declaring the inputs and outputs for this cicle
TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )
OUT_FOR=${FOR/.fastq/.CLEANED.fastq}
OUT_REV=${REV/.fastq/.CLEANED.fastq}

### Using bbduk for filtering reads right portion with quality under 20  and/or  total length under 100 and/or entropy less than 0.01
bbduk.sh in1=$FOLDER_INPUT/$FOR in2=$FOLDER_INPUT/$REV out1=processed_FASTQ/$OUT_FOR out2=processed_FASTQ/$OUT_REV ref=$ADAPTERS_FILE tpe tbo qtrim=r trimq=20 minlen=100 entropy=0.01 outm=FASTQ_check/BBDuk/Discarded/$i.txt 2> FASTQ_check/BBDuk/${i/.fastq.gz/.txt}

# Reducing the space occupied by discarded reads
#grep "^@" FASTQ_check/BBDuk/Discarded/$i.txt -A 1 --no-group-separator | sed "s/@/>/g" > FASTQ_check/BBDuk/Discarded/$i.fasta
#gzip FASTQ_check/BBDuk/Discarded/$i.fasta
rm FASTQ_check/BBDuk/Discarded -r

done




####### EXTRA: QUALITY CONTROL AFTER THE CLEANING #######

#mkdir FASTQ_check/fastqc_processed_reads

#for i in ./processed_FASTQ/*CLEANED_* ; do fastqc $i -o FASTQ_check/fastqc_processed_reads --threads $threads; done

#multiqc FASTQ_check/fastqc_processed_reads -o FASTQ_check/fastqc_processed_reads

#rm FASTQ_check/fastqc_processed_reads/*fastqc.zip


#open FASTQ_check/fastqc_processed_reads/multiqc_report.html
#cat FASTQ_check/fastqc_processed_reads/multiqc_data/multiqc_general_stats.txt | cut -f 1,3,5,7,9 | sed 's/FastQC_mqc-generalstats-fastqc-//g' > FASTQ_check/fastqc_processed_reads/General_info_Sequences_recap_AFTER_THE_CLEANING.txt
#cat FASTQ_check/fastqc_processed_reads/General_info_Sequences_recap_AFTER_THE_CLEANING.txt



# fastqc does not work by remote (?) but I still need this info
for i in ./processed_FASTQ/*gz
do
echo "$i : $(zcat $i | grep "^@" -c)" >> FASTQ_check/remaining_number_of_seqs_after_process.txt
done
sed -e "s#./raw_FASTQ_dir/##g" -e "s/.gz//g" -i FASTQ_check/remaining_number_of_seqs_after_process.txt



zcat processed_FASTQ/mock_AGS_DNAseq_R1_001.CLEANED.fastq.gz | grep "^@" -A 1 --no-group-separator | grep -v "^@" | awk '{print $0, length}' | cut -f 2 -d ' ' | sort | uniq -c > FASTQ_check/Lenghts_distribution_after_BBDUK_R1.txt

zcat processed_FASTQ/mock_AGS_DNAseq_R2_001.CLEANED.fastq.gz | grep "^@" -A 1 --no-group-separator | grep -v "^@" | awk '{print $0, length}' | cut -f 2 -d ' ' | sort | uniq -c > FASTQ_check/Lenghts_distribution_after_BBDUK_R2.txt

