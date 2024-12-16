#!/user/bin/env bash

# Operative System: Ubuntu (Linux)
# required: conda, fastqc, multiqc, bbmap (bbduk.sh) , kraken2.1.3, braken



####### SETTING VARIABLES


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows loads the environment variable "into the script"

### BBDUK ADAPTERS FILE
ADAPTERS_FILE=~/Desktop/programs/bbmap_v39/resources/adapters.fa

### KRAKEN2 DATABASE (removing eukaryotic reads)
k_decont_db="/media/matteo/SSD1/reference/kraken2_enviroment_contam_db/"

### KRAKEN2 DATABASE (classification)
kraken_db='/media/matteo/SSD1/reference/kraken2_microbial_db/'
# inside the same path of kraken database has been built also the braken database (with kmer=35)

### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### THREADS TO USE
threads=31
# nproc # to check the PC availability





####### PREPARING THE DIRECTORIES #######

mkdir processed_FASTQ 
mkdir FASTQ_check
mkdir FASTQ_check/fastqc_raw_reads
mkdir FASTQ_check/BBDuk
mkdir FASTQ_check/BBDuk/Discarded
mkdir FASTQ_check/kraken2_eukar_contam
mkdir kraken_output




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
grep "^@" FASTQ_check/BBDuk/Discarded/$i.txt -A 1 --no-group-separator | sed "s/@/>/g" > FASTQ_check/BBDuk/Discarded/$i.fasta
gzip FASTQ_check/BBDuk/Discarded/$i.fasta
rm FASTQ_check/BBDuk/Discarded/$i.txt
#rm FASTQ_check/BBDuk/Discarded/$i.fasta

done




####### REMOVING EUKAR READS THROUGH KRAKEN2 #######

# the K2 confidence threshold has been chosen comparing the estimated abs human read amount between kraken2 and bowtie2 (the last is much more accurate but it's difficult to handle with huge db!)
# (necessary to avoid problems due to databases imperfections, e.g. see DOI:10.1101/2023.07.28.550993 and DOI:10.1186/s13059-020-02023-1 )

conda activate kraken2


FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep ".CLEANED.fastq" | extract_name )


echo -e "Reads NOT assigned to human, plants or insects ... \n" > FASTQ_check/kraken2_eukar_contam/Not_eukar_percents.txt

mkdir temp_decont


for i in $SAMPLE_LIST
do

echo -e "\n\n\n *** Decontaminating the sample $i with Kraken... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

kraken2 --threads $threads --db $k_decont_db --output temp_decont/${i}_output.txt --report temp_decont/${i}_report.txt --confidence 0.45 --minimum-hit-groups 3 --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV --unclassified-out temp_decont/${i}_#.fastq
# NB: the '#' in the command above is not a comment, is required by kraken for paired reads!

echo -e "in sample $i : $(grep "unclassified" temp_decont/${i}_report.txt | head -n1 | cut -f 1)" >> FASTQ_check/kraken2_eukar_contam/Not_eukar_percents.txt

mv temp_decont/${i}_report.txt FASTQ_check/kraken2_eukar_contam/${i}_report.txt
rm temp_decont/${i}_output.txt

gzip temp_decont/${i}__1.fastq temp_decont/${i}__2.fastq
mv temp_decont/${i}__1.fastq.gz $FOLDER_INPUT/${i}_R1_CLEANED_Decontam.fastq.gz
mv temp_decont/${i}__2.fastq.gz $FOLDER_INPUT/${i}_R2_CLEANED_Decontam.fastq.gz

done

rm -r temp_decont


# Removing the NOT decontam fastq
if [ $(ls $FOLDER_INPUT | grep -c "_Decontam" ) -eq $(ls $FOLDER_INPUT | grep -c ".CLEANED.fastq" ) ]
then
for r in $(ls $FOLDER_INPUT/ | grep -v '_Decontam' ) ; do rm $FOLDER_INPUT/$r ; done
fi


### Generating a summary of contaminants

echo -e "# NB: the abundances of the single organisms do not sum up due to LCA algorithm ( many tax nodes between that level and the domain level) and due to the protozoa which are not listed here" > FASTQ_check/kraken2_eukar_contam/Summary_of_contaminants.tsv
echo -e "Sample\tPercents_on_total\tRead_assigned\tOrganism" >> FASTQ_check/kraken2_eukar_contam/Summary_of_contaminants.tsv

rm -f Summary_temp.tsv
for x in Eukaryota Insecta sapiens Mus Platyhelminthes Nematoda Anellida Chlorophyta Viridiplantae Toxoplasma
do
grep -w $x FASTQ_check/kraken2_eukar_contam/*_report.txt | sed "s#.*/##g" | sed "s/_report.txt//g" | cut -f 1,2,6 | sed "s/   //g" | sed "s/: /:\t/g " | sed "s/Eukaryota/Every_eukaryota/g" >> Summary_temp.tsv
done
cat Summary_temp.tsv | sort >> FASTQ_check/kraken2_eukar_contam/Summary_of_contaminants.tsv
rm Summary_temp.tsv


# conda deactivate





####### STARTING THE ANALYSIS WITH KRAKEN2 #######

# conda activate kraken2

FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "_Decontam" | extract_name )


### Using Kraken2 on each file
for i in $SAMPLE_LIST
do

echo -e "\n\n\n\n *** Processing the sample $i with Kraken... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

kraken2 --threads $threads --db $kraken_db --output kraken_output/${i}_output.txt --report kraken_output/${i}_report.txt --confidence 0.05 --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV

bracken -r 150 -d $kraken_db -i kraken_output/${i}_report.txt -o kraken_output/${i}_bracken_abund.tsv
# -r is the reads length

done


# kraken2-inspect --db $kraken_db --threads 30 | cut -f 4,6 > FASTQ_check/Unclassified_Kraken/Actual_db.tsv


conda deactivate



# searching for unclassified reads ...

# NB : modificare da fastq a fastq.gz e farne un loop

#mkdir FASTQ_check/Unclassified_Kraken
#grep "^U" kraken_output/SRR17089433_output.txt | cut -f 2 > FASTQ_check/unclassified_SRR17089433.txt
#grep -F -f FASTQ_check/Unclassified_Kraken/unclassified_SRR17089433.txt processed_FASTQ/SRR17089433_1.CLEANED.fastq -A1 --no-group-separator | cut -f 1 -d ' ' | sed "s/@/>/g" > FASTQ_check/Unclassified_Kraken/unclassified_SRR17089433.fasta   # NB: do not remove the capital F, otherwise grep will be slow
#rm FASTQ_check/Unclassified_Kraken/unclassified_SRR17089433.txt
#gzip FASTQ_check/Unclassified_Kraken/unclassified_SRR17089433.fasta



# optimizing the required space to store files
rm kraken_output/*_output.txt   # it's a SAM-like (big and no more required for this script)


# from kraken output to phyloseq format
conversion_script="/home/matteo/Desktop/Tool_box/Scripts/Scripts_Shotgun_DNA/Script_from_KrakenBrakenOutput_TO_phyloseq.R"
Rscript $conversion_script




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

