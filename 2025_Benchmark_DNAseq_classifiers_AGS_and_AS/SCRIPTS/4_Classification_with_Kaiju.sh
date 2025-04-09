#!/user/bin/env bash

# required: bash, Kaiju v1.10, Kaiju official database (nt_euk 2023), custom Kaiju db (see methods) and R
# kaiju pre-indexed databases are available at https://bioinformatics-centre.github.io/kaiju/downloads.html

#PATH="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin:$PATH"




####### PREPARING THE ENVIRONMENT #####


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases


### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### THREADS TO USE
threads=33
# nproc # to check the PC availability

### WHERE IS KAIJU
alias kaiju="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju"
alias kaiju2table="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju2table"

### KAIJU DATABASES
Kaiju_DB=/media/matteo/SSD4/reference/kaiju_db_nr_euk_2023-05-10
# Kaiju_DB_no_Euka=/media/matteo/SSD4/reference/kaiju_db_nr_NO_EUKA_2023_05_10
Kaiju_DB_plus=/media/matteo/SSD4/reference/Kaiju_complete_nr_2024


mkdir Kaiju_output
mkdir EXTRA_missclassif_kaiju_mock
mkdir EXTRA_missclassif_kaiju_mock/db_default
mkdir EXTRA_missclassif_kaiju_mock/db_plus    # the complete custom db dedicated folder




####### STARTING THE ANALYSIS WITH KAIJU #######

FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep ".fastq" | extract_name )


### Using Kaiju on each file ...
for i in $SAMPLE_LIST     # main loop
do
TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )



for c in 0.01  0.00001   # nested loop 
# (-e stands for E values, the default is 0.01, the second values comes from https://doi.org/10.3390/d15070835)

do



### Nested-nested loop on pre-index db ...

for m in 11 30 42  # nested-nested loop : testing also higher minimum alignments (default is 11 aa --> 33 nt, hence increasing it at each iteration)
do

echo -e "\n\n\n\n *** Processing the sample $i with Kaiju pre-indexed db, E-value $c and m $m ... *** \n"

kaiju -t $Kaiju_DB/nodes.dmp -f $Kaiju_DB/kaiju_db_*.fmi  -E $c -m $m  -i $FOLDER_INPUT/$FOR -j $FOLDER_INPUT/$REV  -o Kaiju_output/${i}_Evalue${c}_m${m}_output.txt  -z $threads  -x -v
# -x enables SEG algorithm for filtering low-complexity region (suggested by the official tutorial)
# -v enables the addiction of other columns to the output


# EXTRA: checking the mock homo sapiens missclassification ...

# NB: from the official Kaiju GitHub ... "nr_euk --- Like option -s nr and additionally include proteins from fungi and microbial eukaryotes, see taxon list in bin/kaiju-taxonlistEuk.tsv" --> there is NO Homo expected in this db!
cat Kaiju_output/mock_AGS_DNAseq_Evalue${c}_m${m}_output.txt | grep "Homo" | cut -f 2,3 > EXTRA_missclassif_kaiju_mock/db_default/Taxa_classified_from_Homo_reads_conf${c}_m${m}.txt
grep -vw "0" EXTRA_missclassif_kaiju_mock/db_default/Taxa_classified_from_Homo_reads_conf${c}_m${m}.txt > EXTRA_missclassif_kaiju_mock/db_default/Homo_reads_missclassified_as_bacteria_evalue${c}_m${m}.txt
cat EXTRA_missclassif_kaiju_mock/db_default/Homo_reads_missclassified_as_bacteria_evalue${c}_m${m}.txt | grep -Ew "2015173|2015172|55799|93504|29056|51953|51952|5761|5762|1418015|1418016|3764|74632|5820|5854|168256|168256|1437326|1673731|6668|35525|5833|5855|36330|5858|864142" > EXTRA_missclassif_kaiju_mock/db_default/Unexpected_yet_frequent_taxa_classified_conf${c}_m${m}.txt


# Extra check on Diploscapter misclassified as bacteria (not featured in this db)
cat Kaiju_output/mock_AGS_DNAseq_Evalue${c}_m${m}_output.txt | grep "Diploscapter" | grep -vw "0" | cut -f 2,3 > EXTRA_missclassif_kaiju_mock/db_default/Diploscapter_reads_as_bacteriaconf${c}_m${m}.txt


done   # end of the nested nested loop (m, pre-indexed db)




### Nested-nested loop (2) on custom db ...

for m in 11 30 42  # nested-nested loop : testing also higher minimum alignments (default is 11 aa --> 33 nt, hence increasing it at each iteration)
do

echo -e "\n\n\n\n *** Processing the sample $i with Kaiju custom db, E-value $c and m $m ... *** \n"

kaiju -t $Kaiju_DB_plus/nodes.dmp -f $Kaiju_DB_plus/kaiju_db_*.fmi  -E $c -m $m  -i $FOLDER_INPUT/$FOR -j $FOLDER_INPUT/$REV  -o Kaiju_output/${i}_Evalue${c}_m${m}_output_PLUS.txt  -z $threads  -x -v

# Homo is here included ...
cat Kaiju_output/mock_AGS_DNAseq_Evalue${c}_m${m}_output_PLUS.txt | grep "Homo" | cut -f 2,3 > EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Homo_reads_conf${c}_m${m}.txt
cat Kaiju_output/mock_AGS_DNAseq_Evalue${c}_m${m}_output_PLUS.txt  | grep -v "Homo" | cut -f 2,3 | grep -w "9606" > EXTRA_missclassif_kaiju_mock/db_plus/Not_Homo_classified_as_homo_conf$c_m${m}.txt

echo -e "Number of homo reads: $(cat EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Homo_reads_conf${c}_m${m}.txt | wc -l ) \nNumber of TRUE homo reads not recognised: $( grep -v "9606" EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Homo_reads_conf${c}_m${m}.txt | wc -l ) \nNumber of NON homo reads classified as such: $( cat EXTRA_missclassif_kaiju_mock/db_plus/Not_Homo_classified_as_homo_conf$c_m${m}.txt | wc -l ) \n Number of TRUE homo reads properly recognised: $( grep -w "9606" EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Homo_reads_conf${c}_m${m}.txt | wc -l ) " >  EXTRA_missclassif_kaiju_mock/db_plus/Counts_of_homo_false_and_true_positives_conf${c}_m${m}.txt

cat EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Homo_reads_conf${c}_m${m}.txt | grep -Ew "2015173|2015172|55799|93504|29056|51953|51952|5761|5762|1418015|1418016|3764|74632|5820|5854|168256|168256|1437326|1673731|6668|35525|5833|5855|36330|5858|864142" > EXTRA_missclassif_kaiju_mock/db_plus/Unexpected_yet_frequent_taxa_classified_conf${c}_m${m}_FROM_HOMO.txt


# Extra check on that Diploscapter (also included in the plus db)
cat Kaiju_output/mock_AGS_DNAseq_Evalue${c}_m${m}_output_PLUS.txt | grep "Diploscapter" | cut -f 2,3 > EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Diploscapter_reads_conf${c}_m${m}.txt
cat Kaiju_output/mock_AGS_DNAseq_Evalue${c}_m${m}_output_PLUS.txt  | grep -v "Diploscapter" | cut -f 2,3 | grep -Ew "55799|288516|1182519|367193|2018661|2621930|55800" > EXTRA_missclassif_kaiju_mock/db_plus/Not_Diploscapter_classified_as_Diploscapter_conf$c_m${m}.txt

echo -e "Number of Diploscapter reads: $(cat EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Diploscapter_reads_conf${c}_m${m}.txt | wc -l ) \nNumber of NON Diploscapter reads classified as such: $( cat EXTRA_missclassif_kaiju_mock/db_plus/Not_Diploscapter_classified_as_Diploscapter_conf$c_m${m}.txt | wc -l ) \nNumber of TRUE Diploscapter reads properly recognised: $( grep -Ew "55799|288516|1182519|367193|2018661|2621930|55800" EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Diploscapter_reads_conf${c}_m${m}.txt | wc -l ) \nNumber of reads classified as "nematodes" (LCA): $( grep -w "6231" EXTRA_missclassif_kaiju_mock/db_plus/Taxa_classified_from_Diploscapter_reads_conf${c}_m${m}.txt | wc -l ) " >  EXTRA_missclassif_kaiju_mock/db_plus/Counts_of_Diploscapter_false_and_true_positives_conf${c}_m${m}.txt


done   # end of the nested nested loop 2 (m, custom db)



done   # end of the nested loop (E value)


# ### Testing the Kaiju output if the database does not contains also eukaryotes (checking the missclassification of bacteria)
# kaiju -t $Kaiju_DB_no_Euka/nodes.dmp -f $Kaiju_DB_no_Euka/kaiju_db_*.fmi  -E 0.0001 -m 30  -i $FOLDER_INPUT/$FOR -j $FOLDER_INPUT/$REV  -o Kaiju_output/${i}_Evalue_0.001_m22__NO_EUKA_output.txt  -z $threads  -x -v


done   # end of the main loop (samples)


 

### Building a summary of results    # !!! IT WILL BE DONE AFTERWARDS (AFTER THE CONTIGS CLASSIFICATION, SEE SCRIPT 6)
#output_names=$(ls Kaiju_output | grep "output.txt" | sed 's#^#Kaiju_output/#g')  # this is the list of inputs for the row below
#kaiju2table -t $Kaiju_DB/nodes.dmp -n $Kaiju_DB/names.dmp -r genus -l genus -o Kaiju_output/kaiju_classif_summary.tsv $output_names
# inputs written as last argument
# -r choose the depth of taxonomic classification
# -m filters out taxa (of "r" level") with relative abundance below m threshold
# -l chooses which taxon level will be printed in the output
# -p adds the full taxonomic name (from phylum to species), not just the identified one (can't be used with l)

# Few entries are missing because they lack of the genus level in NCBI (but not in SILVA)! Taking these rows from the species level.
#kaiju2table -t $Kaiju_DB/nodes.dmp -n $Kaiju_DB/names.dmp -r species -l species -o Kaiju_output/kaiju_classif_summary_species_level.tsv $output_names


#grep "Moranbacteria bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  > missing_genus_to_add_from_species.txt
# NB: Moranii is Moranbacteria in Kaiju db, differently from Kraken2 nt_core 
#grep "Kaiserbacteria bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
#grep "Nomurabacteria bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
#grep "SH-PL14" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
#grep "Magasanik" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt   # none
#grep "Saccharimonadia bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
#grep "Solirubrobacterales bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt   # 67-14
#grep -e "OPS17|OPS 17" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt # none
#grep -e "37-13" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt # none
#grep -e "Actinobacteria bacterium IMCC26207" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt # none (NB: different codes)
## Transforming the species name (various species) ensuring the manual parsing between SILVA and NCBI databases (which is coded in another script)
#sed "s/Actinobacteria bacterium IMCC26207/IMCC26207/g" -i missing_genus_to_add_from_species.txt
#sed "s/ bacterium.*/ bacterium/g" -i missing_genus_to_add_from_species.txt
#sed "s/uncultured //g" -i missing_genus_to_add_from_species.txt
#cat Kaiju_output/kaiju_classif_summary.tsv  missing_genus_to_add_from_species.txt > temp.tsv   # both of the files (original and updated) in an unique file
#mv temp.tsv Kaiju_output/kaiju_classif_summary.tsv   # overwriting the original one
#rm missing_genus_to_add_from_species.txt

# NB: now different rows of the same species are in the genus output !




#### EXTRA: Translating the tax ID codes of the extra check files of mock miss classifications...

sed -i "s/51953/...missclassif as... Elaeis guineensis/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/2015173/...missclassif as...Ooceraea biroi/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/2015172/...missclassif as...Ooceraea/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/55799/...missclassif as...Diploscapter/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/55800/...missclassif as...Diploscapter sp./g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/93504/...missclassif as...Ostrinia furnacalis/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/29056/...missclassif as...Ostrinia/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/51953/...missclassif as...Elaeis guineensis/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/51952/..missclassif as...Elaeis/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5761/...missclassif as...Naegleria/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5762/...missclassif as...Naegleria gruberi/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1418015/...missclassif as...Pseudochloris/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1418016/...missclassif as...Pseudochloris/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/3764/...missclassif as...Rosa/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/74632/...missclassif as...Rosa gallica/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5820/...missclassif as...Plasmodium/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5854/...missclassif as...Plasmodium reichenowi/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5833/...missclassif as...Plasmodium falciparum/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/5855/...missclassif as...Plasmodium vivax/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/36330/...missclassif as...Plasmodium ovale/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/864142/...missclassif as...Plasmodium ovale wallikeri/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/\b5858\b/...missclassif as...Plasmodium malarie/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/168256/...missclassif as...Zoothamnium/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1437326/...missclassif as...Zoothamnium intermedium/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/1673731/...missclassif as...Suigetsumonas/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/\b6668\b/...missclassif as...Daphnia/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*
sed -i "s/35525/...missclassif as...Daphnia magna/g" EXTRA_missclassif_kaiju_mock/*/Unexpected_yet_frequent_taxa_classified_conf*



# rm Kaiju_output/*_output.txt    # cleaning
### These files are still required ... cleaning in the 6th script, after the parsing of Kaiju's files (including contigs classif)




############# EXTRA: CHECKING THE EFFICIENCY OF VIRUS DEDICATED DB ##################

#### rvdb v29
# Where_db="/media/matteo/SSD4/reference"
# wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_rvdb_2024-12-20.tgz -O $Where_db/rvdb_v29
# mkdir $Where_db/kaiju_db_rvdb_v29
# tar -zxvf $Where_db/rvdb_v29 -C $Where_db/kaiju_db_rvdb_v29

#kaiju_viral_db=$Where_db/Kaiju_viral_refseq_2024"
#FORW_INPUT=processed_FASTQ/mock_AGS_DNAseq_R1_001.CLEANED.fastq.gz
#REV_INPUT=processed_FASTQ/mock_AGS_DNAseq_R2_001.CLEANED.fastq.gz
#kaiju -t $kaiju_viral_db/nodes.dmp -f $kaiju_viral_db/kaiju_db_*.fmi  -E 0.01 -m 22  -i $FORW_INPUT -j $REV_INPUT -o Kaiju_output/Check_rvdb29_output.txt  -z $threads  -x -v
#grep "phage" Kaiju_output/Check_rvdb29_output.txt | grep "^C" | wc -l > Number_of_phageT4_reads_classified.txt
#rm Kaiju_output/Check_rvdb29_output.txt 
# *** ANY T4 SEQUENCE WAS CLASSIFIED AT ALL! ***



#### viral refseq
Where_db="/media/matteo/SSD4/reference"
#wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_viruses_2024-08-15.tgz -O $Where_db/Kaiju_viral_refseq
#mkdir $Where_db/Kaiju_viral_refseq_2024
#tar -zxvf $Where_db/Kaiju_viral_refseq -C $Where_db/Kaiju_viral_refseq_2024

kaiju_viral_db=$Where_db/Kaiju_viral_refseq_2024
FORW_INPUT=processed_FASTQ/mock_AGS_DNAseq_R1_001.CLEANED.fastq.gz
REV_INPUT=processed_FASTQ/mock_AGS_DNAseq_R2_001.CLEANED.fastq.gz
for m in 25 40
do
kaiju -t $kaiju_viral_db/nodes.dmp -f $kaiju_viral_db/kaiju_db_*.fmi  -E 0.01 -m $m  -i $FORW_INPUT -j $REV_INPUT -o Kaiju_output/Check_kaiju_viral_output.txt  -z 40  -x -v
cut -f 1,2,3 Kaiju_output/Check_kaiju_viral_output.txt | grep "phage" | grep "^C" > EXTRA_missclassif_kaiju_mock/Reads_T4_classified_with_RefSeq_viral_only.txt
echo -e "Reads of T4 phage correctly classified from Kaiju on virus RefSeq: $(grep -Ew "2681598|10665|10663" EXTRA_missclassif_kaiju_mock/Reads_T4_classified_with_RefSeq_viral_only.txt | wc -l ) \nReads (T4 and not) missclassfied: $(grep -v "phage" Kaiju_output/Check_kaiju_viral_output.txt | grep -Ev "2681598|10665|10663" | wc -l) " > EXTRA_missclassif_kaiju_mock/Counts_T4_RefSeq_viral_m${m}.txt
rm  Kaiju_output/Check_kaiju_viral_output.txt  EXTRA_missclassif_kaiju_mock/Reads_T4_classified_with_RefSeq_viral_only.txt
done



### can blastN detect the seq?
#zcat mock_AGS_DNAseq_R1_001.CLEANED.fastq.gz | grep "NC_000866.4_Enterobacteria_phage_T4,_1229_62" -A 1
#@NC_000866.4_Enterobacteria_phage_T4,_1229_62/1
#ATCCATATACTTAAATGCTTCTGTCAACTGCGATTTAAGGCATTCGCAAATTGAAAGAGAATTTTTTGTATTAGGTTTAAAGCTAATTCGTGCAGTTCTATTATTTTCTTTTAACGGTCTTACTTCCATCTGAAGAGTATAACCATCA
# BlastN --> "somewhat dissimilar sequences"
# perfect matches with "Escherichia phages"   --> yes!

