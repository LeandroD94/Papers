##### START ####
perl scripts/chrom_collect.pl lavori_nostri/A-Lavoro_FAS_Biopsie_Feci_Saliva_CRC/bowtie2_SAM_output.txt lavori_nostri/A-Lavoro_FAS_Biopsie_Feci_Saliva_CRC/FAS_counts_OTU_biopsia_contaminata_con_NA.csv FAS_biopsia > analyzed/FAS_biopsia.txt
perl scripts/chrom_collect.pl lavori_nostri/A-Lavoro_FAS_Biopsie_Feci_Saliva_CRC/bowtie2_SAM_output.txt lavori_nostri/A-Lavoro_FAS_Biopsie_Feci_Saliva_CRC/counts_OTU_Saliva_FAS.csv > analyzed/FAS_saliva.txt
perl scripts/chrom_collect.pl lavori_nostri/A-Lavoro_FAS_Biopsie_Feci_Saliva_CRC/bowtie2_SAM_output.txt lavori_nostri/A-Lavoro_FAS_Biopsie_Feci_Saliva_CRC/counts_OTU_Stool_FAS.csv > analyzed/FAS_stool.txt

perl scripts/chrom_collect.pl lavori_nostri/B-Valvole_aortiche_CAVD/bowtie2_SAM_output.txt lavori_nostri/B-Valvole_aortiche_CAVD/counts_otu_pre_decontaminaz.csv > analyzed/CAVD.txt






#### END #####



#this transforms fasta into two columns records
cd reofftargethostdurante16sseq
cat Every_ASV_Sequence.fasta | perl -ne 's/.\n$/ /;print' | perl -ne 's/>/\n>/g;print' > Every_ASV_Sequence.fasta.tab



cd Collezione_dati_HOST_off_target_2_The_revenge

cd "A-Lavoro_FAS_Biopsie_Feci_Saliva_CRC"
perl contam_count.pl eucar_contaminant_ASV.fasta FAS_counts_OTU_biopsia_contaminata_con_NA.csv > read_counts_human_biopsia_plus.txt
perl contam_count.pl eucar_contaminant_ASV.fasta counts_OTU_Saliva_FAS.csv > read_counts_human_saliva_plus.txt
perl contam_count.pl eucar_contaminant_ASV.fasta counts_OTU_Stool_FAS.csv > read_counts_human_stool_plus.txt
cd ..
cd "B-Valvole_aortiche  CAVD"
perl contam_count.pl eucar_contaminant_ASV.fasta counts_otu_pre_decontaminaz.csv > read_counts_human_CAVD_plus.txt
cd ..
cd "C-ACC_Surrenali_sane_e_Tumorali"
perl contam_count.pl eucar_contaminant_ASV.fasta ASV_raw_abundances_phyla_without_decontam.csv > read_counts_human_ACC_plus.txt
cd ..
cd "D-I_leggendari_batteri_nel_sangue"
perl contam_count.pl eucar_contaminant_ASV.fasta



grep chr bowtie2_SAM_output.txt | cut -f 1,3 > chroms.txt


perl kmer_screen.pl eucar_contaminant_ASV.fasta 14 100 > kmer14_screen_out.txt
perl kmer_screen.pl eucar_contaminant_ASV.fasta 15 100 > kmer15_screen_out.txt
perl kmer_screen.pl eucar_contaminant_ASV.fasta 16 100 > kmer16_screen_out.txt
perl kmer_screen.pl eucar_contaminant_ASV.fasta 17 100 > kmer17_screen_out.txt

head -n6 read_counts_human_biopsia.txt
d763d514fa6a4d6b3810535b169c3305	21301
9e62fb049a7e60bb1f0241449ac244d9	11768
a520b40f23556b7546b3df416e302d05	6050
c0586f98e28a03ace4a7d4442f15f4ad	3539
38d8a32ebd53ed2a523d45d314b4a582	1874
6d3b195e7894967193f59ca85df78b0c	1802

head -n6 read_counts_human_biopsia.txt | cut -f1
d763d514fa6a4d6b3810535b169c3305
9e62fb049a7e60bb1f0241449ac244d9
a520b40f23556b7546b3df416e302d05
c0586f98e28a03ace4a7d4442f15f4ad
38d8a32ebd53ed2a523d45d314b4a582
6d3b195e7894967193f59ca85df78b0c

#FASgrep.pl -l reads_most_abundant.txt eucar_contaminant_ASV.fasta

FASgrep.pl d763d514fa6a4d6b3810535b169c3305 eucar_contaminant_ASV.fasta > reads_most_abundant.fasta
FASgrep.pl 9e62fb049a7e60bb1f0241449ac244d9 eucar_contaminant_ASV.fasta >> reads_most_abundant.fasta
FASgrep.pl a520b40f23556b7546b3df416e302d05 eucar_contaminant_ASV.fasta >> reads_most_abundant.fasta
FASgrep.pl c0586f98e28a03ace4a7d4442f15f4ad eucar_contaminant_ASV.fasta >> reads_most_abundant.fasta
FASgrep.pl 38d8a32ebd53ed2a523d45d314b4a582 eucar_contaminant_ASV.fasta >> reads_most_abundant.fasta
FASgrep.pl 6d3b195e7894967193f59ca85df78b0c eucar_contaminant_ASV.fasta >> reads_most_abundant.fasta

perl kmer_screen.pl reads_most_abundant.fasta 14 100 > kmer14_screen_most_abundant_out.txt
perl kmer_screen.pl reads_most_abundant.fasta 15 100 > kmer15_screen_most_abundant_out.txt
perl kmer_screen.pl reads_most_abundant.fasta 16 100 > kmer16_screen_most_abundant_out.txt
perl kmer_screen.pl reads_most_abundant.fasta 17 100 > kmer17_screen_most_abundant_out.txt

#final selection: the following maps the top 4 most abundant reads
head -n3 kmer17_screen_most_abundant_out.txt | tail -n1
TGATAAACCTTTAGCAA	4	55.27	4	>c0586f98e28a03ace4a7d4442f15f4ad|>a520b40f23556b7546b3df416e302d05|>9e62fb049a7e60bb1f0241449ac244d9|>d763d514fa6a4d6b3810535b169c3305


