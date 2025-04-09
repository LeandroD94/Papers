#!/user/bin/env bash

# Operative System: Ubuntu (Linux)
# required: conda, MEGAHIT 1.2.9, MetaBat 2:2.17 ,  kMetaShot v1
# official kMerShot databases downloaded from   https://benlangmead.github.io/aws-indexes/k2   on feb 2025



####### SETTING VARIABLES

### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows loads the environment variable "into the script"

### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### THREADS TO USE
threads=42
# nproc # to check the PC availability

### DATABASE KMETASHOT
kmetashot_db='/media/matteo/SSD4/reference/kMetaShot_reference.h5'




############# CREATING AND CLASSIFYING THE MAGs ###############


FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "CLEANED" | extract_name )

mkdir MAGs_output/


for i in $SAMPLE_LIST      # nested loop (on each sample)
do

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )


### Using MEGAHIT on each file to create contigs

conda activate MEGAHIT

# NB: MEGAHIT uses by default the >90% of the PC available resources

echo -e "\n\n\n\n *** Processing the sample $i with MEGAHIT (default settings)... *** \n"
megahit -1 $FOLDER_INPUT/$FOR -2 $FOLDER_INPUT/$REV -o MAGs_output/${i}_contigs_default 

echo -e "\n\n\n\n *** Processing the sample $i with MEGAHIT (preset meta-large)... *** \n"
megahit -1 $FOLDER_INPUT/$FOR -2 $FOLDER_INPUT/$REV --presets meta-large -o MAGs_output/${i}_contigs_metalarge
# from the help: meta-large is designed for large & complex metagenomes, like soil 
# ... in practice, the min kmer is increased compared to the default, the highest is decreased and the "k step" is set to 10

echo -e "\n\n\n\n *** Processing the sample $i with MEGAHIT (custom klist from MetaShot paper)... *** \n"
megahit -1 $FOLDER_INPUT/$FOR -2 $FOLDER_INPUT/$REV --k-list 35,57,79,99 -o MAGs_output/${i}_contigs_custom 
# the k min employed by the MetaShot paper is much higher!

rm MAGs_output/${i}_contigs_*/intermediate_contigs* -r



### ASSEMBLING CONTIGS IN MAGs WITH METABAT2
#   !!!  same custom parameters as kMetaShot !!! which are --maxP 99 --minS 98 -m 1,500
# -m [ --minContig ] (default=2500)    Minimum size of a contig for binning
# --maxP (default=95)                  Percentage of 'good' contigs considered for binning decided by connection among contigs. The greater, the more sensitive.
# --minS (default=60)                  Minimum score of a edge for binning (should be between 1 and 99). The greater, the more specific.
# NB: the "sensitivity" is required for simple communities, according to the metabat1 manual, see   https://gensoft.pasteur.fr/docs/MetaBAT/2.15/ the help, they 
# -t [ --numThreads ] (default=0, meaning that all of the cores will be used)  

# NB: same conda environment as MEGAHIT


echo -e "\n\n\n\n *** Assembling the sample $i contigs in MAGs with MetaBat2 ... *** \n"

mkdir MAGs_output/${i}_contigs_default/MAGs
metabat -i MAGs_output/${i}_contigs_default/final.contigs.fa -o MAGs_output/${i}_contigs_default/MAGs/MAG_default --maxP 99 --minS 98 -m 1500 --seed 1994

mkdir MAGs_output/${i}_contigs_metalarge/MAGs
metabat -i MAGs_output/${i}_contigs_metalarge/final.contigs.fa -o MAGs_output/${i}_contigs_metalarge/MAGs/MAGs_metalarge --maxP 99 --minS 98 -m 1500 --seed 1994

mkdir MAGs_output/${i}_contigs_custom/MAGs
metabat -i MAGs_output/${i}_contigs_custom/final.contigs.fa -o MAGs_output/${i}_contigs_custom/MAGs/MAGs_custom --maxP 99 --minS 98 -m 1500 --seed 1994




### CLASSIFYING THE MAGs WITH kMETASHOT #######

conda deactivate  # changing env
conda activate kmetashot

threads=45  # <- !
# !!!!! NB: if the number of bins or headers is lower than the number of processors then
# an error like "Shape of passed values is (x, y), indices imply (x, z)" will be obtained !!!


for c in 0.2 0.4   # confidence levels (nested loop)
do

echo -e "\n\n\n\n *** Classifying the sample $i contigs with kMETASHOT (conf $c)... *** \n"
kMetaShot_classifier_NV.py -b MAGs_output/${i}_contigs_default/final.contigs.fa  -r $kmetashot_db   -p $threads -o MAGs_output/${i}_contigs_classif_default_conf${c}  -a $c > temp.log
kMetaShot_classifier_NV.py -b MAGs_output/${i}_contigs_metalarge/final.contigs.fa  -r $kmetashot_db   -p $threads -o MAGs_output/${i}_contigs_classif_metalarge_conf${c}  -a $c > temp.log
kMetaShot_classifier_NV.py -b MAGs_output/${i}_contigs_custom/final.contigs.fa  -r $kmetashot_db   -p $threads -o MAGs_output/${i}_contigs_classif_custom_conf${c}  -a $c > temp.log

rm MAGs_output/${i}_contigs_classif*/bins -r
rm MAGs_output/${i}_contigs_classif*/tmp -r

echo -e "\n\n\n\n *** Classifying the sample $i MAGs with kMETASHOT (conf $c)... *** \n"
kMetaShot_classifier_NV.py -b MAGs_output/${i}_contigs_default/MAGs  -r $kmetashot_db   -p $threads -o MAGs_output/${i}_MAGs_classif_default_conf${c}  -a $c > temp.log
kMetaShot_classifier_NV.py -b MAGs_output/${i}_contigs_metalarge/MAGs  -r $kmetashot_db   -p $threads -o MAGs_output/${i}_MAGs_classif_metalarge_conf${c}  -a $c > temp.log
kMetaShot_classifier_NV.py -b MAGs_output/${i}_contigs_custom/MAGs  -r $kmetashot_db   -p $threads -o MAGs_output/${i}_MAGs_classif_custom_conf${c}  -a $c > temp.log

rm MAGs_output/${i}_MAGs_classif*/bins -r

rm temp.log   # it was required just because kMetaShot is VERY verbose

# cut  -f 5,6,12 MAGs_output/${i}_contigs__classif_conf${c}/kMetaShot_classification_resume_a2r* -d , | grep -Ew "2015173|2015172|55799|93504|29056|51953|51952|5761|5762|1418015|1418016|3764|74632|5820|5854|168256|168256|1437326|1673731|6668|35525|5833|5855|36330|5858" > MAGs_output/${i}_contigs__classif_conf${c}/Unexpected_yet_frequent_taxa_missclassified_conf$c.txt


done # end of the nested loop (conf)

done # end of the main loop (sample)

conda deactivate




### EXTRA: CLASSIFYING THE CONTINGS WITH KAIJU AND KRAKEN2 #######


threads=60
kraken_db='/media/matteo/SSD4/reference/kraken2_core_nt_db_v20241228/'
alias kaiju="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju"
alias kaiju2table="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju2table"
#Kaiju_DB=/media/matteo/SSD4/reference/kaiju_db_nr_euk_2023-05-10
Kaiju_DB=/media/matteo/SSD4/reference/Kaiju_complete_nr_2024   # using the "full" complete database


### Using Kraken2

conda activate kraken2

echo -e "\n\n\n\n *** Processing the contigs (metalarge) with Kraken (Conf $c) ... *** \n"

c=0.15 # kraken confidence
mkdir kraken_output/output_contigs_conf$c
kraken2 --threads $threads --db $kraken_db --output kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_output.txt --report kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_report.txt --confidence $c  MAGs_output/mock_AGS_DNAseq_contigs_metalarge/final.contigs.fa
bracken -r 150 -d $kraken_db -i kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_report.txt  -l G  -o kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Genus.tsv  #NB: genus level
bracken -r 150 -d $kraken_db -i kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_report.txt  -l S  -o kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  #NB: species level

# Few entries does not have any "genus" name according to NCBI, hence it is not featured in the genus report... but they are abundant in the reactor according to the 16S dataset --> adding them manually from the species level table
grep "Moraniibacteriota bacterium" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  > missing_genus_to_add_from_species.txt
grep "Kaiserbacteria bacterium" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Nomurabacteria bacterium" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "SH-PL14" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Magasanik" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt   # none
grep "Saccharimonadia bacterium" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Solirubrobacterales bacterium" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt   # 67-14
grep -e "OPS17|OPS 17" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none
grep -e "37-13" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none
grep -e "Actinobacteria bacterium IMCC26207" kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none (NB: different codes)
cat kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Genus.tsv  missing_genus_to_add_from_species.txt > temp.tsv   # both of the files (original and updated) in an unique file
mv temp.tsv   kraken_output/output_contigs_conf$c/mock_AGS_DNAseq_bracken_abund_Genus.tsv  # overwriting the original one
rm missing_genus_to_add_from_species.txt

conda deactivate


### Using Kaiju...
echo -e "\n\n\n\n *** Processing the contigs (metalarge) with Kaiju ... *** \n"
kaiju -t $Kaiju_DB/nodes.dmp -f $Kaiju_DB/kaiju_db_*.fmi  -E 0.0001 -m 42  -i MAGs_output/mock_AGS_DNAseq_contigs_metalarge/final.contigs.fa  -o Kaiju_output/mock_AGS_DNAseq_contigs_Evalue0.0001_m42_output_PLUS.txt  -z $threads  -x -v




### END
# saving space after the script ...
gzip MAGs_output/${i}_contigs_*/final.contigs.fa
gzip MAGs_output/${i}_contigs_*/MAGs -r



