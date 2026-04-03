#!/user/bin/env bash

# Operative System: Ubuntu (Linux)
# required: conda, bbmap (bbduk.sh), Kaiju v1.10



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
threads=15
# nproc # to check the PC availability

### KAIJU VARIABLES
#PATH="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin:$PATH"
alias kaiju="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju"
alias kaiju2table="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju2table"
Kaiju_DB=/media/matteo/SSD4/reference/Kaiju_complete_nr_2024  # custom complete database
Kaiju_to_Phylo_script=/home/matteo/Desktop/Tool_box/Scripts/Scripts_DNAseq_microb/Script_from_KaijuSummary_TO_phyloseq.R




####### PREPARING THE DIRECTORIES #######

mkdir processed_FASTQ 
mkdir FASTQ_check
# mkdir FASTQ_check/fastqc_raw_reads
mkdir FASTQ_check/BBDuk
mkdir FASTQ_check/BBDuk/Discarded
mkdir Kaiju_output




####### INITIAL QUALITY CONTROL #######

#for i in ./raw_FASTQ_dir/*gz ; do fastqc $i -o FASTQ_check/fastqc_raw_reads/ --threads $threads ;  echo -e "\n"  ; done

#multiqc FASTQ_check/fastqc_raw_reads/ -o FASTQ_check/fastqc_raw_reads

#open FASTQ_check/fastqc_raw_reads/multiqc_report.html
#cat FASTQ_check/fastqc_raw_reads/multiqc_data/multiqc_general_stats.txt | cut -f 1,3,5,7,9 | sed 's/FastQC_mqc-generalstats-fastqc-//g' > #FASTQ_check/fastqc_raw_reads/General_info_Sequences_recap.txt

#cat FASTQ_check/fastqc_raw_reads/General_info_Sequences_recap.txt

#rm FASTQ_check/fastqc_raw_reads/*fastqc.zip



# fastqc does not work by remote (?) but I still need this info
for i in ./raw_FASTQ_dir/*gz
do 
echo "$i : $(zcat $i | grep "^@" -c)" >> FASTQ_check/original_number_of_seqs.txt
done
sed -e "s#./raw_FASTQ_dir/##g" -e "s/.gz//g" -i FASTQ_check/original_number_of_seqs.txt




####### REMOVING THE BAD SEQUENCES WITH BBDUK #######

FOLDER_INPUT="raw_FASTQ_dir"
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



####### STARTING THE ANALYSIS WITH KAIJU #######

FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep ".CLEANED.fastq" | extract_name )


### Using Kaiju on each file ...
for i in $SAMPLE_LIST     # main loop
do
TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

m=30      # min length threshold  (NB: few bioprojects are about 100 bp!)

echo -e "\n\n\n\n *** Processing the sample $i with Kaiju custom db,  m= $m ... *** \n"

kaiju -t $Kaiju_DB/nodes.dmp -f $Kaiju_DB/kaiju_db_plus.fmi  -m $m  -i $FOLDER_INPUT/$FOR -j $FOLDER_INPUT/$REV  -o Kaiju_output/${i}_m${m}_output.txt  -z $threads  -x -v
# -x enables SEG algorithm for filtering low-complexity region (suggested by the official tutorial)
# -v enables the addiction of other columns to the output

done  


### Building a summary of results
output_names=$(ls Kaiju_output | grep "output.txt" | sed 's#^#Kaiju_output/#g')  # this is the list of inputs for the row below
kaiju2table -p -r species -t $Kaiju_DB/nodes.dmp -n $Kaiju_DB/names.dmp -e -o Kaiju_output/kaiju_classif_summary.tsv $output_names
# inputs written as last argument
# -r choose the depth of taxonomic classification
# -m filters out taxa (of "r" level") with relative abundance below m threshold
# -l chooses which taxon level will be printed in the output
# -p adds the full taxonomic name (from phylum to species), not just the identified one (can't be used with l)
# -e is required to maintain the viruses which are not classifiable down to the taxonomic level specified with "r"

# kaiju2table -t $Kaiju_DB/nodes.dmp -n $Kaiju_DB/names.dmp -r genus -l genus -o Kaiju_output/kaiju_classif_summary_genus_level.tsv $output_names

gzip Kaiju_output/kaiju_classif_summary.tsv
#rm Kaiju_output/*output*


# building the phyloseq object
Rscript $Kaiju_to_Phylo_script




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
sed -e "s#./processed_FASTQ/##g" -e "s/.gz//g" -i FASTQ_check/remaining_number_of_seqs_after_process.txt


#zcat processed_FASTQ/mock_AGS_DNAseq_R1_001.CLEANED.fastq.gz | grep "^@" -A 1 --no-group-separator | grep -v "^@" | awk '{print $0, length}' | cut -f 2 -d ' ' | sort | uniq -c > FASTQ_check/Lenghts_distribution_after_BBDUK_R1.txt
#zcat processed_FASTQ/mock_AGS_DNAseq_R2_001.CLEANED.fastq.gz | grep "^@" -A 1 --no-group-separator | grep -v "^@" | awk '{print $0, length}' | cut -f 2 -d ' ' | sort | uniq -c > FASTQ_check/Lenghts_distribution_after_BBDUK_R2.txt

