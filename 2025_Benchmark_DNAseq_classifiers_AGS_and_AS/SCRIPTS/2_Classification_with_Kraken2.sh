#!/user/bin/env bash

# Operative System: Ubuntu (Linux)
# required: conda, kraken2.1.3, braken
# official kraken2 databases downloaded from   https://benlangmead.github.io/aws-indexes/k2   on feb 2025



####### SETTING VARIABLES

### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows loads the environment variable "into the script"

### KRAKEN2 DATABASE
kraken_db='/media/matteo/SSD4/reference/kraken2_core_nt_db_v20241228/'

### KRAKEN2 - SILVA DATABASE
kraken_SILVA='/media/matteo/SSD4/reference/kraken2_on_SILVA138_db/'

### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### THREADS TO USE
threads=62
# nproc # to check the PC availability




####### STARTING THE ANALYSIS WITH KRAKEN2  #######

conda activate kraken2

FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "CLEANED" | extract_name )

mkdir EXTRA_missclassif_kraken_mock
mkdir kraken_output/


for c in 0.05 0.15 0.3 0.45 0.65 0.85 0.99  # main loop (kraken confidence levels)
do

mkdir kraken_output/output_conf$c/



### Using Kraken2 on each file
for i in $SAMPLE_LIST      # nested loop (on each sample)
do

echo -e "\n\n\n\n *** Processing the sample $i with Kraken (Conf $c) ... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

# processing using nt_core database...
kraken2 --threads $threads --db $kraken_db --output kraken_output/output_conf$c/${i}_output.txt --report kraken_output/output_conf$c/${i}_report.txt --confidence $c --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV

bracken -r 150 -d $kraken_db -i kraken_output/output_conf$c/${i}_report.txt  -l G  -o kraken_output/output_conf$c/${i}_bracken_abund_Genus.tsv  #NB: genus level
bracken -r 150 -d $kraken_db -i kraken_output/output_conf$c/${i}_report.txt  -l S  -o kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  #NB: species level
# -r is the reads length


# Few entries does not have any "genus" name according to NCBI, hence it is not featured in the genus report... but they are abundant in the reactor according to the 16S dataset --> adding them manually from the species level table
grep "Moraniibacteriota bacterium" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  > missing_genus_to_add_from_species.txt
grep "Kaiserbacteria bacterium" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Nomurabacteria bacterium" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "SH-PL14" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Magasanik" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt   # none
grep "Saccharimonadia bacterium" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Solirubrobacterales bacterium" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt   # 67-14
grep -e "OPS17|OPS 17" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none
grep -e "37-13" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none
grep -e "Actinobacteria bacterium IMCC26207" kraken_output/output_conf$c/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none (NB: different codes)
cat kraken_output/output_conf$c/${i}_bracken_abund_Genus.tsv  missing_genus_to_add_from_species.txt > temp.tsv   # both of the files (original and updated) in an unique file
mv temp.tsv kraken_output/output_conf$c/${i}_bracken_abund_Genus.tsv  # overwriting the original one
rm missing_genus_to_add_from_species.txt


# checking the mock homo sapiens missclassification ...
if $( echo $i | grep -q "mock_" )  # if the sample is the mock then ...
then

cat kraken_output/output_conf$c/mock_AGS_DNAseq_output.txt | grep "Homo" | cut -f 2,3 > EXTRA_missclassif_kraken_mock/Taxa_classified_from_Homo_reads_conf$c.txt
cat kraken_output/output_conf$c/mock_AGS_DNAseq_output.txt | grep -v "Homo" | cut -f 2,3 | grep -w "9606" > EXTRA_missclassif_kraken_mock/Not_Homo_classified_as_homo_conf$c.txt
echo -e "Number of homo reads: $(cat EXTRA_missclassif_kraken_mock/Taxa_classified_from_Homo_reads_conf$c.txt | wc -l ) \nNumber of TRUE homo reads not recognised: $( grep -v "9606" EXTRA_missclassif_kraken_mock/Taxa_classified_from_Homo_reads_conf$c.txt | wc -l ) \nNumber of NON homo reads classified as such: $( cat EXTRA_missclassif_kraken_mock/Not_Homo_classified_as_homo_conf$c.txt | wc -l )" >  EXTRA_missclassif_kraken_mock/Counts_of_homo_false_and_true_positives_conf$c.txt

cat kraken_output/output_conf$c/mock_AGS_DNAseq_output.txt | grep "phage" | cut -f 2,3 > EXTRA_missclassif_kraken_mock/Taxa_classified_from_T4_reads_conf$c.txt
echo -e "Number of total phage T4 reads in the mock: $(cat EXTRA_missclassif_kraken_mock/Taxa_classified_from_T4_reads_conf$c.txt | wc -l ) \nNumber of TRUE T4 reads not recognised: $( grep -Ev "2681598|10665|10663" EXTRA_missclassif_kraken_mock/Taxa_classified_from_T4_reads_conf$c.txt | wc -l ) " >  EXTRA_missclassif_kraken_mock/Counts_of_T4_false_and_true_positives_conf$c.txt

# further question: does any of these known read classifies as common eukaryotes found in our past dataset?
cat kraken_output/output_conf$c/mock_AGS_DNAseq_output.txt | cut -f 2,3 | grep -Ew "2015173|2015172|55799|93504|29056|51953|51952|5761|5762|1418015|1418016|3764|74632|5820|5854|168256|168256|1437326|1673731|6668|35525|5833|5855|36330|5858" > EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf$c.txt
# 2015173 is Ooceraea biroi, 2015172 is Oocearea (genus level), 55799 is Diploscapter (a nematode, genus level), 93504 is Ostrinia furnacalis (an insect), 29056 is Ostrinia (genus level), 51953 is Elaeis guineensis (a palm),
# 51952 is Elaeis, 5761 is Naegleria (the brain eating amoeba genus), 5762 is Naegleria gruberi, 1418015 and 1418016 are Pseudochloris and its ref species (algae), 3764 is the genus of the roses (flower), 74632 is Rosa gallica (ref species), 
# 5820 is Plasmodium and 5854 is Plasmodium reichenowi, 5833 is P. falciparum, 5855 is P. vivax, 36330 P.ovale, 5858 P.malarie, 
# 168256 is Zoothamnium and 1437326 is Zoothamnium intermedium, 1673731 is Suigetsumonas, 6668 is Daphnia and 35525 is Daphnia magna

cat kraken_output/output_conf$c/mock_AGS_DNAseq_output.txt | grep -v "Plasmodium" | cut -f 2,3 | grep -Ew "5820|5833|5854|5850|5855|5825|5861" > EXTRA_missclassif_kraken_mock/Not_Plasmodium_classified_as_Plasmodium_conf$c.txt

fi  # end of the missclassification check



### again, but this time kraken is run on SILVA138
#echo -e "\n\n *** Again, but with Kraken-SILVA ... *** \n"
#kraken2 --threads $threads --db $kraken_SILVA --output kraken_output/output_conf$c/${i}_SILVA_output.txt --report kraken_output/output_conf$c/${i}_SILVA_report.txt --confidence $c --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV
#bracken -r 150 -d $kraken_SILVA -i kraken_output/output_conf$c/${i}_SILVA_report.txt  -l G  -o kraken_output/output_conf$c/${i}_SILVA_bracken_abund_Genus.tsv  #NB: genus level



# optimizing the required space to store files
rm kraken_output/output_conf$c/*_output.txt   # it's a SAM-like (big and no more required from now on)


done   # nested loop (sample)
done   # main loop (confidence level)


# kraken2-inspect --db $kraken_db --threads 30 | cut -f 4,6 > FASTQ_check/Unclassified_Kraken/Actual_db.tsv


# conda deactivate



# from kraken output to phyloseq format
#conversion_script="/home/matteo/Desktop/Tool_box/Scripts/Scripts_Shotgun_DNA/Script_from_KrakenBrakenOutput_TO_phyloseq.R"
#Rscript $conversion_script




#### EXTRA: Translating the tax ID codes of the extra check files of mock miss classifications...

sed -i "s/51953/...missclassif as... Elaeis guineensis/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/2015173/...missclassif as...Ooceraea biroi/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/2015172/...missclassif as...Ooceraea/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/55799/...missclassif as...Diploscapter/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/55800/...missclassif as...Diploscapter sp./g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/93504/...missclassif as...Ostrinia furnacalis/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/29056/...missclassif as...Ostrinia/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/51953/...missclassif as...Elaeis guineensis/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/51952/..missclassif as...Elaeis/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5761/...missclassif as...Naegleria/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5762/...missclassif as...Naegleria gruberi/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1418015/...missclassif as...Pseudochloris/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1418016/...missclassif as...Pseudochloris/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/3764/...missclassif as...Rosa/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/74632/...missclassif as...Rosa gallica/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5820/...missclassif as...Plasmodium/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5854/...missclassif as...Plasmodium reichenowi/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5833/...missclassif as...Plasmodium falciparum/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5855/...missclassif as...Plasmodium vivax/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/36330/...missclassif as...Plasmodium ovale/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5858/...missclassif as...Plasmodium malarie/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/168256/...missclassif as...Zoothamnium/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1437326/...missclassif as...Zoothamnium intermedium/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1673731/...missclassif as...Suigetsumonas/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/6668/...missclassif as...Daphnia/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/35525/...missclassif as...Daphnia magna/g" EXTRA_missclassif_kraken_mock/Unexpected_yet_frequent_taxa_classified_conf*




#### EXTRA: Translating the tax ID codes of the extra check files of mock miss classifications...

cat EXTRA_missclassif_kraken_mock/Taxa_classified_from_Homo_reads_conf0.15.txt | grep -v "9606" | cut -f 2 > EXTRA_missclassif_kraken_mock/Taxa_list_misclassified_from_Homo_reads_conf0.15.txt
cat EXTRA_missclassif_kraken_mock/Taxa_classified_from_Homo_reads_conf0.45.txt | grep -v "9606" | cut -f 2 > EXTRA_missclassif_kraken_mock/Taxa_list_misclassified_from_Homo_reads_conf0.45.txt
cat EXTRA_missclassif_kraken_mock/Taxa_classified_from_Homo_reads_conf0.99.txt | grep -v "9606" | cut -f 2 > EXTRA_missclassif_kraken_mock/Taxa_list_misclassified_from_Homo_reads_conf0.99.txt

# these files will be used later (script n8)



########### EXTRA: TESTING THE CONFIDENCE REQUIRED BY THE CUSTOM EUKARYOTE DB ############


### KRAKEN2 EUKAR CUSTOM DATABASE (1)
kraken_db_euka='/media/matteo/SSD4/reference/kraken2_enviroment_contam_db/'
# database built with refseqs of hexapoda, anellides, chlorophyta, plants, h sapiens, mus musculus, Univec_core


### KRAKEN2 EUKAR CUSTOM DATABASE (2)
#kraken_db_euka_all='/media/matteo/SSD4/reference/OLD_kraken2_enviroment_contam_db/'
# if the db includes also plathelmintes and nematodes (="euka ALL")...


# conda activate kraken2   # re-activating (as the this eukaryote extra check is optional)

mkdir EXTRA_missclassif_kraken_mock/WHAT_IF_EUKAR_DB
OUT_FOLDER="EXTRA_missclassif_kraken_mock/WHAT_IF_EUKAR_DB"
FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "mock" | grep "CLEANED" | extract_name )


for c in 0.45 0.65 0.85 0.99 # main loop (kraken confidence levels)
do


### Using Kraken2 on each file
for i in $SAMPLE_LIST      # nested loop (on each sample)
do

echo -e "\n\n\n\n *** Processing $i with custom EUKAR Kraken (Conf $c) ... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )


### processing using the database 1...
kraken2 --threads $threads --db $kraken_db_euka --output $OUT_FOLDER/${i}_euka_conf${c}_output.txt --report $OUT_FOLDER/${i}_euka_conf${c}_report.txt --confidence $c --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV

cat $OUT_FOLDER/${i}_euka_conf${c}_output.txt | grep "Homo" | cut -f 2,3 > $OUT_FOLDER/Taxa_classified_from_Homo_reads_conf$c.txt
cat $OUT_FOLDER/${i}_euka_conf${c}_output.txt | grep -v "Homo" | cut -f 2,3 | grep -w "9606" > $OUT_FOLDER/Not_Homo_classified_as_homo_conf$c.txt
echo -e "Number of homo reads: $(cat $OUT_FOLDER/Taxa_classified_from_Homo_reads_conf$c.txt | wc -l ) \nNumber of true homo reads not recognised: $( grep -v "9606" $OUT_FOLDER/Taxa_classified_from_Homo_reads_conf$c.txt | wc -l ) \nNumber of NON homo reads classified as such: $( cat $OUT_FOLDER/Not_Homo_classified_as_homo_conf$c.txt | wc -l )" >  $OUT_FOLDER/Counts_of_homo_false_and_true_positives_conf$c.txt

cat $OUT_FOLDER/${i}_euka_conf${c}_output.txt | cut -f 2,3 | grep -w "10090" > $OUT_FOLDER/Not_Mouse_classified_as_Mouse_conf$c.txt

cat $OUT_FOLDER/${i}_euka_conf${c}_output.txt | cut -f 2,3 | grep -Ew "2015173|2015172|55799|93504|29056|51953|51952|5761|5762|1418015|1418016|3764|74632|5820|5854|168256|168256|1437326|1673731|6668|35525|4650|94328|5833|5855|36330|5858" > $OUT_FOLDER/Unexpected_yet_frequent_taxa_classified_conf$c.txt


### processing using the database 2...
#kraken2 --threads $threads --db $kraken_db_euka_all --output $OUT_FOLDER/${i}__eukaALSOnemat_conf${c}_output.txt --report $OUT_FOLDER/${i}__eukaALSOnemat_conf${c}_report.txt --confidence $c --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV

#cat $OUT_FOLDER/${i}__eukaALSOnemat_conf${c}_output.txt | grep "Homo" | cut -f 2,3 > $OUT_FOLDER/Taxa_classified_from_Homo_reads_eukaALL_conf$c.txt
#cat $OUT_FOLDER/${i}__eukaALSOnemat_conf${c}_output.txt | grep -v "Homo" | cut -f 2,3 | grep -w "9606" > $OUT_FOLDER/Not_Homo_classified_as_homo_eukaALL_conf$c.txt
#echo -e "Number of homo reads: $(cat $OUT_FOLDER/Taxa_classified_from_Homo_reads_eukaALL_conf$c.txt | wc -l ) \nNumber of true homo reads not recognised: $( grep -v "9606" $OUT_FOLDER/Taxa_classified_from_Homo_reads_eukaALL_conf$c.txt | wc -l ) \nNumber of NON homo reads classified as such: $( cat $OUT_FOLDER/Not_Homo_classified_as_homo_eukaALL_conf$c.txt | wc -l )" >  $OUT_FOLDER/Counts_of_homo_false_and_true_eukaALLdb_conf$c.txt

#cat $OUT_FOLDER/${i}__eukaALSOnemat_conf${c}_output.txt | cut -f 2,3 | grep -w "10090" > $OUT_FOLDER/Not_Mouse_classified_as_Mouse_eukaALL_conf$c.txt

#cat $OUT_FOLDER/${i}__eukaALSOnemat_conf${c}_output.txt | cut -f 2,3 | grep -Ew "2015173|2015172|55799|93504|29056|51953|51952|5761|5762|1418015|1418016|3764|74632|5820|5854|168256|168256|1437326|1673731|6668|35525|4650|94328|5833|5855|36330|5858" > $OUT_FOLDER/Unexpected_yet_frequent_taxa_class_eukaALLdb_conf$c.txt


# optimizing the required space to store files
rm $OUT_FOLDER/*_output.txt   # it's a SAM-like (big and no more required from now on)



done   # nested loop (sample)
done   # main loop (confidence level)



### Translating the tax ID codes of the extra check files of mock miss classifications...

sed -i "s/2015173/...missclassif as...Ooceraea biroi/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/2015172/...missclassif as...Ooceraea/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/55799/...missclassif as...Diploscapter/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/55800/...missclassif as...Diploscapter sp./g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/93504/...missclassif as...Ostrinia furnacalis/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/29056/...missclassif as...Ostrinia/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/51953/...missclassif as...Elaeis guineensis/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/51952/..missclassif as...Elaeis/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/5761/...missclassif as...Naegleria/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/5762/...missclassif as...Naegleria gruberi/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/1418015/...missclassif as...Pseudochloris/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/1418016/...missclassif as...Pseudochloris/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/3764/...missclassif as...Rosa/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/74632/...missclassif as...Rosa gallica/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/5820/...missclassif as...Plasmodium/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/5854/...missclassif as...Plasmodium reichenowi/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/5833/...missclassif as...Plasmodium falciparum/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/5855/...missclassif as...Plasmodium vivax/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/36330/...missclassif as...Plasmodium ovale/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/5858/...missclassif as...Plasmodium malarie/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/168256/...missclassif as...Zoothamnium/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/1437326/...missclassif as...Zoothamnium intermedium/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/1673731/...missclassif as...Suigetsumonas/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/6668/...missclassif as...Daphnia/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/35525/...missclassif as...Daphnia magna/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/4650/...missclassif as...Zingiber/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*
sed -i "s/94328/...missclassif as...Zingiber/g" $OUT_FOLDER/Unexpected_yet_frequent_taxa*

conda deactivate

