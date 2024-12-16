# Operative System: Ubuntu (Linux)
# required: conda, bbmap (bbduk.sh) , kraken2.1.3, sortmerna 4.3.6, datasets 16.4.3,  R 4.3
# optional: fastqc, multiqc, bracken, salmon 1.10

# script version 13/08/2024




####### SETTING VARIABLES


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows loads the environment variable "into the script"
alias diamond=~/Desktop/programs/diamond

### BBDUK ADAPTERS FILE
ADAPTERS_FILE=/home/matteo/Desktop/programs/bbmap_v39/resources/adapters.fa

### KRAKEN2 DATABASES
k_decont_db='/media/matteo/SSD1/reference/kraken2_enviroment_contam_db/'
kraken_db='/media/matteo/SSD1/reference/kraken2_microbial_db/'   # optional
# inside the same path of kraken micr db has been built also the bracken db (kmer=35)

### SORTMERNA DATABASES
database_rRNA="/media/matteo/SSD1/reference/sortmerna_db/smr_v4.3_default_db.fasta"
# path to the fasta to index
database_ncRNA="/media/matteo/SSD1/reference/prokaryotes_ncRNA/nc_microb_db_scfbio_240124_only_longer_19nt.fasta"  # sequences >19nt are required
database_path="/media/matteo/SSD1/reference/sortmerna_db/"  # where are the already indexed ones

### SALMON ncRNA DATABASE
ncRNA_db_Salmon="/media/matteo/SSD1/reference/prokaryotes_ncRNA/Salmon_ncRNAmicrob_noRRNA_k27"   # optional

### DIAMOND DATABASE
diamond_db="/media/matteo/SSD1/reference/UniRef_Diamond"
# in the same folder is stored also the references file used from the R parsing script

### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"    # e.g. "sample_R1_001.fastq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### THREADS TO USE
threads=32
# nproc # to check the PC availability




####### PREPARING THE DIRECTORIES #######

mkdir processed_FASTQ 
mkdir FASTQ_check
mkdir FASTQ_check/fastqc_raw_reads
mkdir FASTQ_check/BBDuk
mkdir FASTQ_check/BBDuk/Discarded
mkdir FASTQ_check/kraken2_eukar_contam
mkdir FASTQ_check/SortMeRNA
mkdir FASTQ_check/Salmon_ncRNAmapped_no_rRNA
mkdir FASTQ_check/Diamond_mapping_percentages
mkdir kraken_output
mkdir Diamond_output




####### INITIAL QUALITY CONTROL #######

# for i in ./raw_FASTQ_dir/*gz ; do fastqc $i -o FASTQ_check/fastqc_raw_reads/ --threads $threads;  echo -e "\n"  ; done

# multiqc FASTQ_check/fastqc_raw_reads/ -o FASTQ_check/fastqc_raw_reads

##open FASTQ_check/fastqc_raw_reads/multiqc_report.html
#cat FASTQ_check/fastqc_raw_reads/multiqc_data/multiqc_general_stats.txt | cut -f 1,3,4,5,7,9 | sed 's/FastQC_mqc-generalstats-fastqc-//g' > FASTQ_check/fastqc_raw_reads/General_info_Sequences_recap.txt

#cat FASTQ_check/fastqc_raw_reads/General_info_Sequences_recap.txt

#rm FASTQ_check/fastqc_raw_reads/*fastqc.zip




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
gzip FASTQ_check/BBDuk/Discarded/$i.fasta -f   # -f is "force" (overwrites)
rm FASTQ_check/BBDuk/Discarded/$i.txt
# rm FASTQ_check/BBDuk/Discarded/$i.fasta

# rm -r FASTQ_check/BBDuk/Discarded/  # delete this folder if it is not required (save a lot of GB)

done




####### REMOVING EUKAR READS THROUGH KRAKEN2 #######

# the K2 confidence threshold has been chosen comparing the estimated abs human read amount between kraken2 and bowtie2 (the last is much more accurate but it's difficult to handle with huge db!)


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

conda deactivate




####### CHECKING WITH SALMON THE ncRNA ABUNDANCES (EXCLUDING rRNA) #######
# NB: this is an optional check, not required to get the main results 

conda activate salmon


mkdir Salmon_temp

FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "_Decontam" | extract_name )

for i in $SAMPLE_LIST
do

echo -e "\n\n\n *** Quantifing ncRNA in sample $i with Salmon ... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

salmon quant -p $threads --gcBias --posBias -i $ncRNA_db_Salmon -l A -1 $FOLDER_INPUT/$FOR -2 $FOLDER_INPUT/$REV -o Salmon_temp/$i 2> temp.log

grep -v "0.000" Salmon_temp/$i/quant.sf | cut -f 1,2,5 | sed "s/\..*//g" | sed "s/ncRNA_//g" > FASTQ_check/Salmon_ncRNAmapped_no_rRNA/ncRNA_found_in_$i.tsv
# NB: only the rows with at least one count are saved to save space

done

rm temp.log

conda deactivate


grep -E "percent_mapped" Salmon_temp/*/aux_info/meta_info.json | sed "s#Salmon_temp/##g" | sed "s#/aux.*:##g" > FASTQ_check/Salmon_ncRNAmapped_no_rRNA/report_mapped.txt


rm -r Salmon_temp




####### SEARCHING MICROBES THROUGH KRAKEN2 #######   
# NB: this is an optional check, not required to get the main results

conda activate kraken2


FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "_Decontam" | extract_name )


for i in $SAMPLE_LIST
do

echo -e "\n\n\n *** Processing the sample $i with Kraken... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

kraken2 --threads $threads --db $kraken_db --confidence 0.1 --output kraken_output/${i}_output.txt --report kraken_output/${i}_report.txt --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV

bracken -r 150 -d $kraken_db -i kraken_output/${i}_report.txt -o kraken_output/${i}_braken_abund.tsv

done


conda deactivate


# optimizing the required space to store files
rm kraken_output/*_output.txt   # it's a SAM-like (big and not required)


grep "U    0" kraken_output/*report.txt | cut -f 1,2  > FASTQ_check/Reads_UNclassified_as_microbial_by_Kraken2.tsv




####### CLEANING SAMPLES FROM EVERY NCRNA WITH SORTMERNA ####### 

mkdir SortMeRNA   # it's a temp dir required by the program


FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "_Decontam.fastq" | extract_name )


echo "Number of reads identified as ncRNA:" > FASTQ_check/SortMeRNA/Percentuals_of_ncRNA.txt
for i in $SAMPLE_LIST    ### beginning of the loop
do

echo -e "\n\n\n *** Removing every ncRNA from sample $i with SortMeRNA... *** \n"

### declaring the inputs and outputs for this cicle
TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

# using both rRNA db and every other ncRNA db  
sortmerna  --ref $database_rRNA --ref $database_ncRNA --reads $FOLDER_INPUT/$FOR --reads $FOLDER_INPUT/$REV  --workdir SortMeRNA  -fastx -other -paired_in --idx-dir $database_path  --index 0 --out2 --threads $threads > log_sortmerna.txt
# NB: 'index 0' prevents generation of databases, remove it if new ones are needed

rm log_sortmerna.txt
rm -rf SortMeRNA/kvdb/ SortMeRNA/readb  # have to be deleted manually for each run!

echo "$(grep "passing" SortMeRNA/out/aligned.log | sed "s/.*=//")"  # to display it
echo "$i = $(grep "passing" SortMeRNA/out/aligned.log | sed "s/.*=//")" >> FASTQ_check/SortMeRNA/Percentuals_of_ncRNA.txt

### exporting the outputs
mv SortMeRNA/out/aligned.log FASTQ_check/SortMeRNA/aligned_$i.log
mv SortMeRNA/out/other_fwd.fq.gz $FOLDER_INPUT/${FOR/.fastq/_OUTncRNA.fastq}
mv SortMeRNA/out/other_rev.fq.gz $FOLDER_INPUT/${REV/.fastq/_OUTncRNA.fastq}

done

rm -r SortMeRNA



### Removing the fastq with rRNA
if [ $(ls $FOLDER_INPUT | grep -c "_OUTncRNA" ) -eq $(ls $FOLDER_INPUT | grep -c "_Decontam.fastq" ) ]
then
for r in $(ls $FOLDER_INPUT/ | grep -v '_OUTncRNA' ) ; do rm $FOLDER_INPUT/$r ; done
fi




####### MAPPING TO UNIREF USING DIAMOND AND THEN QUANTIFY THE RESULT ####### 

# to use paired information through the so called 'congruent strategy'
diamond_as_paired="/home/matteo/Desktop/Tool_box/Scripts/Scripts_RNAseq_microbial_community/Diamond_as_paired.R"

# to count mapping among every sample
parsing_diamond="/home/matteo/Desktop/Tool_box/Scripts/Scripts_RNAseq_microbial_community/Parsing_Diamond_mappings.R"


FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "_OUTncRNA" | extract_name )

for i in $SAMPLE_LIST
do

echo -e "\n\n\n *** Translating and mapping cDNA of sample $i ... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

echo -e "\n ... Forward ... \n"
diamond blastx -d $diamond_db/Diamond_Uniref100.dmnd -q $FOLDER_INPUT/$FOR -o Diamond_output/matches_1.tsv --evalue 0.00001 --outfmt 6 qseqid sseqid slen evalue --unal 1 -b 6 -c 1 --quiet
echo -e "\n ... Reverse ... \n"
diamond blastx -d $diamond_db/Diamond_Uniref100.dmnd -q $FOLDER_INPUT/$REV -o Diamond_output/matches_2.tsv --evalue 0.00001 --outfmt 6 qseqid sseqid slen evalue --unal 1 -b 6 -c 1 --quiet
# NB: reduce -b and increase -c (or use the default values) if the RAM does not support the operation with those settings

echo -e "\n Evaluating both mappings... \n"
Rscript $diamond_as_paired  Diamond_output  $i

mv Diamond_output/Matching_percent_$i.tsv  FASTQ_check/Diamond_mapping_percentages/
rm Diamond_output/matches_*  # Diamond concatenates if same name!

done


# re-ordering the matching percents files in an unique file
cat FASTQ_check/Diamond_mapping_percentages/Matching_percent_* > FASTQ_check/Diamond_mapping_percentages/Matching_percentS.tsv
sort FASTQ_check/Diamond_mapping_percentages/Matching_percentS.tsv -r | uniq > FASTQ_check/Diamond_matching_percents.tsv
rm FASTQ_check/Diamond_mapping_percentages -R


# parsing the outputs
Rscript $parsing_diamond   Diamond_output  'microb_only'  $diamond_db


# compressing the simil-BAM to save space
tar -zcf Diamond_mappings.tar.gz Diamond_output
rm Diamond_output/Diamond_*

mv  *_counts_from_Diamond.tsv  Diamond_output
mv  IDs_infos.tsv  Diamond_output
mv Diamond_mappings.tar.gz Diamond_output/DiamondMappings.tar.gz

cp FASTQ_check/Diamond_matching_percents.tsv Diamond_output/   # this file is important for further analyses --> putting it also in the main results folder




####### QUICK PRELIMINARY ANALYSIS OF THE RESULTS #######

quick_analysis_script="/home/matteo/Desktop/Tool_box/Scripts/Scripts_RNAseq_microbial_community/Prelim_analysis_of_RNAseq_microbial_Diamond.R"

Rscript $quick_analysis_script  Diamond_output




########## EXTRA: COUNTING EUKARYOTES READS ##########

eukar_collecting_script="/home/matteo/Desktop/Tool_box/Scripts/Scripts_RNAseq_microbial_community/Gathering_eukar_counts_from_Kraken2_and_Diamond.R"

Rscript $eukar_collecting_script




####### EXTRA: QUALITY CONTROL AFTER THE CLEANING #######

#mkdir FASTQ_check/fastqc_processed_reads

#for i in ./processed_FASTQ/*_OUTncRNA* ; do fastqc $i -o FASTQ_check/fastqc_processed_reads --threads $threads; done

#multiqc FASTQ_check/fastqc_processed_reads -o FASTQ_check/fastqc_processed_reads

#rm FASTQ_check/fastqc_processed_reads/*fastqc.zip


## open FASTQ_check/fastqc_processed_reads/multiqc_report.html
#cat FASTQ_check/fastqc_processed_reads/multiqc_data/multiqc_general_stats.txt | cut -f 1,3,4,5,7,9 | sed 's/FastQC_mqc-generalstats-fastqc-//g' > FASTQ_check/fastqc_processed_reads/General_info_Sequences_recap_AFTER_THE_CLEANING.txt
#cat FASTQ_check/fastqc_processed_reads/General_info_Sequences_recap_AFTER_THE_CLEANING.txt

