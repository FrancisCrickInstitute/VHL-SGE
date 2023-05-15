#20230408_VHL_final_pipeline.sh
#Date: 20230408
#This is the set of scripts used to process sequencing data to variant counts
#Meant to be run on each of 6 sets of fastq.gz files, each set in a different folder corresponding to different sets of experiments.
#Not coded to be run start-to-end as a single command, but interactively.

##### MODULES REQUIRED
#  SeqPrep
#  needleall, from EMBOSS package, from GCC/5.4.0-2.26
#  Python 2.7.5
#
#

##### CUSTOM SCRIPTS REQUIRED
# ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py    			#will produce seqprep scripts to run on samples
# ~/home/users/findlag/bin/run_remove_n_bases.py 					#will produce remove_n_bases scripts to run for all samples
# ~/home/users/findlag/bin/remove_n_bases.py 						#will remove reads from fastq files if N bases present
# ~/home/users/findlag/bin/run_cDNA_to_gDNA_SGE_pipeline.py 		#will produce scripts to run below code for all RNA samples
# ~/home/users/findlag/bin/cDNA_to_gDNA_SGE_pipeline.py 			#will convert cDNA reads to gDNA reads to match sequence context for comparison
# ~/home/users/findlag/bin/run_needle_to_sam_sbatch_pipeline.py 	#will produce needleall scripts to run on samples for alignment
# ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py 	#will perform pairwise analysis of cigar frequencies in 2 samples
# ~/home/users/findlag/bin/20220826_sam_to_edits_v9s_pipeline.py 	#will count hdr events and SNV edits for SGE experiments
# ~/home/users/findlag/bin/220407_annotate_variants_pipeline.py 	#will annotate variants with CADD data, ClinVar data (where available)
# ~/home/users/bucklem/_Scripts/_VHL/Analysis_pipelines/pad_to_gDNA_VHLx1c-hdrl.py 	#will pad target fastq file with adapter ref sequence on end

##### input files requires (in addition to fastqs)
# /camp/lab/findlayg/home/shared/projects/SGE/VHL/fasta
# /camp/lab/findlayg/home/shared/projects/SGE/VHL/VHL_editing_data.txt
# /camp/lab/findlayg/home/shared/projects/SGE/VHL/VHL_cDNA_data.txt
# /camp/lab/findlayg/home/shared/projects/SGE/VHL/cadd
# /camp/lab/findlayg/home/shared/projects/SGE/BRCA1/amplicon_coords_hg19.txt
# /camp/lab/findlayg/home/shared/projects/SGE/VHL/ClinVar_VHL_210618_Single_nucleotide_copy.txt

##### a folder containing necessary scripts and input files (apart from fastq's):  /camp/lab/findlayg/home/shared/projects/SGE/VHL/VHL_sh_pipeline_20240411


#Defining parameters for all samples across 
adapter1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" #this is the nextera reverse adapter, as seen in R1
adapter2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"  #this is the nextera forward adapter, as seen in R2
seq_type="N" #for nextseq
ref_folder="/camp/lab/findlayg/home/shared/projects/SGE/VHL/fasta"
editing_info_file="/camp/lab/findlayg/home/shared/projects/SGE/VHL/VHL_editing_data.txt"
cDNA_info_file="/camp/lab/findlayg/home/shared/projects/SGE/VHL/VHL_cDNA_data.txt"
reads_threshold=".000002"
alignment_score_threshold="300"
orientation="+"
#Consensus coding sequence for gene, should match what's used in ClinVar
CCDS="CCDS2597.1"
amplicon_coords="/camp/lab/findlayg/home/shared/projects/SGE/BRCA1/amplicon_coords_hg19.txt"
#directory with cadd files for each amplicon, named as {amplicon}.cadd e.g. BRCA1x18.cadd or X2.cadd
cadd_dir="/camp/lab/findlayg/home/shared/projects/SGE/VHL/cadd"
#Must download a file from clinvar with all variants selected (e.g. at least 1 star on this date, and then convert it's format to tab-delimited txt file and reference below)
cadd_version="cadd.1.6"
clinvar_file="/camp/lab/findlayg/home/shared/projects/SGE/VHL/ClinVar_VHL_210618_Single_nucleotide_copy.txt"
min_indel_freq=".00002"

#directory in which to run pipeline, will be where all fastq files go and pipeline outputs
my_dir="/camp/lab/findlayg/home/shared/projects/SGE/VHL/20230408_VHL_all"

mkdir /camp/lab/findlayg/home/shared/projects/SGE/VHL/20230408_VHL_all
mkdir /camp/lab/findlayg/home/shared/projects/SGE/VHL/20230408_VHL_all/fastq
mkdir /camp/lab/findlayg/home/shared/projects/SGE/VHL/20230408_VHL_all/pipeline_output/
mkdir /camp/lab/findlayg/home/shared/projects/SGE/VHL/20230408_VHL_all/pipeline_output/sam
mkdir /camp/lab/findlayg/home/shared/projects/SGE/VHL/20230408_VHL_all/pipeline_output/cigar_counts
mkdir /camp/lab/findlayg/home/shared/projects/SGE/VHL/20230408_VHL_all/pipeline_output/editing_data
mkdir $my_dir/pipeline_output/variants


# LOCATION OF ALL FASTQ FILES --->>>>  $my_dir/fastq    [to contain 6 folders, one for each date]

#####FOR EACH ITERATION corresponding to a set of experiments (date), MUST CHANGE:
#dates in folder paths
#cigar_comparisons to match experiments
#amplicon_list to match experiments
#exp_groupings to match experiments



##### 20210722

cigar_comparisons="VHLx3a_HDRL+VHLx3a_neg,VHLx3a_HDRL+VHLx3ar1_D6,VHLx3ar1_D6+VHLx3ar1_D13,VHLx3ar1_D6+VHLx3ar1_D20,VHLx3ar1_D13+VHLx3ar1_D20,VHLx3a_HDRL+VHLx3a_neg,VHLx3a_HDRL+VHLx3ar2_D6,VHLx3ar2_D6+VHLx3ar2_D13,VHLx3ar2_D6+VHLx3ar2_D20,VHLx3ar2_D13+VHLx3ar2_D20,VHLx3b_HDRL+VHLx3b_neg,VHLx3b_HDRL+VHLx3br1_D6,VHLx3br1_D6+VHLx3br1_D13,VHLx3br1_D6+VHLx3br1_D20,VHLx3br1_D13+VHLx3br1_D20,VHLx3b_HDRL+VHLx3b_neg,VHLx3b_HDRL+VHLx3br2_D6,VHLx3br2_D6+VHLx3br2_D13,VHLx3br2_D6+VHLx3br2_D20,VHLx3br2_D13+VHLx3br2_D20"
amplicon_list="VHLx3a,VHLx3b" #splits on commas
#exp name + 8 sample entries per grouping to be uniform across all runs (some experiments have 5 samples, some 7) must be specified, order must be:  ( experiment + pre + post + HDRLlib + neg + rna + post2 + post3 + rna2 , ... )
exp_groupings="VHLx3arL41+VHLx3ar1_D6+VHLx3ar1_D20+VHLx3a_HDRL+VHLx3a_neg+VHLx3ar1_D13+VHLx3ar1_D13+VHLx3ar1_D13+VHLx3ar1_D13,VHLx3arL42+VHLx3ar2_D6+VHLx3ar2_D20+VHLx3a_HDRL+VHLx3a_neg+VHLx3ar2_D13+VHLx3ar2_D13+VHLx3ar2_D13+VHLx3ar2_D13,VHLx3brL41+VHLx3br1_D6+VHLx3br1_D20+VHLx3b_HDRL+VHLx3b_neg+VHLx3br1_D13+VHLx3br1_D13+VHLx3br1_D13+VHLx3br1_D13,VHLx3brL42+VHLx3br2_D6+VHLx3br2_D20+VHLx3b_HDRL+VHLx3b_neg+VHLx3br2_D13+VHLx3br2_D13+VHLx3br2_D13+VHLx3br2_D13"


mkdir $my_dir/pipeline_output/sam/20210722
cd $my_dir/fastq/20210722

mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns
module load GCC/5.4.0-2.26
module load EMBOSS

mkdir sam
python ~/home/users/findlag/bin/run_needle_to_sam_sbatch_pipeline.py $ref_folder
sbatch run_needle_to_sam_sbatch.sh 

mv sam/* $my_dir/pipeline_output/sam/20210722
cd $my_dir/pipeline_output/sam/20210722

mkdir cigar_counts
python ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -20 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt

mkdir variant_counts
python /camp/lab/findlayg/home/users/findlag/bin/20220826_sam_to_edits_v9s_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold
cd variant_counts
mkdir final
python ~/home/users/findlag/bin/220407_annotate_variants_pipeline.py $amplicon_list $ref_folder $editing_info_file $orientation $CCDS $amplicon_coords $cadd_dir $cadd_version $clinvar_file 
cd final
head -10 *summary.txt >> combined_editing_data.txt

mkdir $my_dir/pipeline_output/editing_data/20210722
mv combined_editing_data.txt $my_dir/pipeline_output/editing_data/20210722
mv *summary.txt $my_dir/pipeline_output/editing_data/20210722
mkdir $my_dir/pipeline_output/variants/20210722
mv *.txt $my_dir/pipeline_output/variants/20210722
mkdir $my_dir/pipeline_output/cigar_counts/20210722
mv ../../cigar_counts/*.txt $my_dir/pipeline_output/cigar_counts/20210722
cd ../..
rmdir cigar_counts
rmdir variant_counts/final
rm variant_counts/*.txt
rmdir variant_counts


##### 20210729

cigar_comparisons="VHLx1a_HDRL+VHLx1a_neg,VHLx1a_HDRL+VHLx1ar1_D6,VHLx1ar1_D6+VHLx1ar1_D13,VHLx1ar1_D6+VHLx1ar1_D20,VHLx1ar1_D13+VHLx1ar1_D20,VHLx1a_HDRL+VHLx1a_neg,VHLx1a_HDRL+VHLx1ar2_D6,VHLx1ar2_D6+VHLx1ar2_D13,VHLx1ar2_D6+VHLx1ar2_D20,VHLx1ar2_D13+VHLx1ar2_D20,VHLx1b_HDRL+VHLx1b_neg,VHLx1b_HDRL+VHLx1br1_D6,VHLx1br1_D6+VHLx1br1_D13,VHLx1br1_D6+VHLx1br1_D20,VHLx1br1_D13+VHLx1br1_D20,VHLx1b_HDRL+VHLx1b_neg,VHLx1b_HDRL+VHLx1br2_D6,VHLx1br2_D6+VHLx1br2_D13,VHLx1br2_D6+VHLx1br2_D20,VHLx1br2_D13+VHLx1br2_D20,VHLx1c_HDRL+VHLx1c_neg,VHLx1c_HDRL+VHLx1cr1_D6,VHLx1cr1_D6+VHLx1cr1_D13,VHLx1cr1_D6+VHLx1cr1_D20,VHLx1cr1_D13+VHLx1cr1_D20,VHLx1c_HDRL+VHLx1c_neg,VHLx1c_HDRL+VHLx1cr2_D6,VHLx1cr2_D6+VHLx1cr2_D13,VHLx1cr2_D6+VHLx1cr2_D20,VHLx1cr2_D13+VHLx1cr2_D20,VHLx2_HDRL+VHLx2_neg,VHLx2_HDRL+VHLx2r1_D6,VHLx2r1_D6+VHLx2r1_D13,VHLx2r1_D6+VHLx2r1_D20,VHLx2r1_D13+VHLx2r1_D20,VHLx2_HDRL+VHLx2_neg,VHLx2_HDRL+VHLx2r2_D6,VHLx2r2_D6+VHLx2r2_D13,VHLx2r2_D6+VHLx2r2_D20,VHLx2r2_D13+VHLx2r2_D20,VHLx1p_HDRL+VHLx1p_neg,VHLx1p_HDRL+VHLx1pr1_D6,VHLx1pr1_D6+VHLx1pr1_D13,VHLx1pr1_D6+VHLx1pr1_D20,VHLx1pr1_D13+VHLx1pr1_D20,VHLx1p_HDRL+VHLx1p_neg,VHLx1p_HDRL+VHLx1pr2_D6,VHLx1pr2_D6+VHLx1pr2_D13,VHLx1pr2_D6+VHLx1pr2_D20,VHLx1pr2_D13+VHLx1pr2_D20"
amplicon_list="VHLx1a,VHLx1b,VHLx1c,VHLx1p,VHLx2" #splits on commas
#exp name + 8 sample entries per grouping to be uniform across all runs (some experiments have 5 samples, some 7) must be specified, order must be:  ( experiment + pre + post + HDRLlib + neg + rna + post2 + post3 + rna2 , ... )
exp_groupings="VHLx1arL41+VHLx1ar1_D6+VHLx1ar1_D20+VHLx1a_HDRL+VHLx1a_neg+VHLx1ar1_D13+VHLx1ar1_D13+VHLx1ar1_D13+VHLx1ar1_D13,VHLx1arL42+VHLx1ar2_D6+VHLx1ar2_D20+VHLx1a_HDRL+VHLx1a_neg+VHLx1ar2_D13+VHLx1ar2_D13+VHLx1ar2_D13+VHLx1ar2_D13,VHLx1brL41+VHLx1br1_D6+VHLx1br1_D20+VHLx1b_HDRL+VHLx1b_neg+VHLx1br1_D13+VHLx1br1_D13+VHLx1br1_D13+VHLx1br1_D13,VHLx1brL42+VHLx1br2_D6+VHLx1br2_D20+VHLx1b_HDRL+VHLx1b_neg+VHLx1br2_D13+VHLx1br2_D13+VHLx1br2_D13+VHLx1br2_D13,VHLx1crL41+VHLx1cr1_D6+VHLx1cr1_D20+VHLx1c_HDRL+VHLx1c_neg+VHLx1cr1_D13+VHLx1cr1_D13+VHLx1cr1_D13+VHLx1cr1_D13,VHLx1crL42+VHLx1cr2_D6+VHLx1cr2_D20+VHLx1c_HDRL+VHLx1c_neg+VHLx1cr2_D13+VHLx1cr2_D13+VHLx1cr2_D13+VHLx1cr2_D13,VHLx1prL41+VHLx1pr1_D6+VHLx1pr1_D20+VHLx1p_HDRL+VHLx1p_neg+VHLx1pr1_D13+VHLx1pr1_D13+VHLx1pr1_D13+VHLx1pr1_D13,VHLx1prL42+VHLx1pr2_D6+VHLx1pr2_D20+VHLx1p_HDRL+VHLx1p_neg+VHLx1pr2_D13+VHLx1pr2_D13+VHLx1pr2_D13+VHLx1pr2_D13,VHLx2rL41+VHLx2r1_D6+VHLx2r1_D20+VHLx2_HDRL+VHLx2_neg+VHLx2r1_D13+VHLx2r1_D13+VHLx2r1_D13+VHLx2r1_D13,VHLx2rL42+VHLx2r2_D6+VHLx2r2_D20+VHLx2_HDRL+VHLx2_neg+VHLx2r2_D13+VHLx2r2_D13+VHLx2r2_D13+VHLx2r2_D13"

mkdir $my_dir/pipeline_output/sam/20210729
cd $my_dir/fastq/20210729

mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns
module load GCC/5.4.0-2.26
module load EMBOSS

mkdir sam
python ~/home/users/findlag/bin/run_needle_to_sam_sbatch_pipeline.py $ref_folder
sbatch run_needle_to_sam_sbatch.sh 

mv sam/* $my_dir/pipeline_output/sam/20210729
cd $my_dir/pipeline_output/sam/20210729

mkdir cigar_counts
python ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -20 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt

mkdir variant_counts
python /camp/lab/findlayg/home/users/findlag/bin/20220826_sam_to_edits_v9s_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold
cd variant_counts
mkdir final
python ~/home/users/findlag/bin/220407_annotate_variants_pipeline.py $amplicon_list $ref_folder $editing_info_file $orientation $CCDS $amplicon_coords $cadd_dir $cadd_version $clinvar_file 
cd final
head -10 *summary.txt >> combined_editing_data.txt

mkdir $my_dir/pipeline_output/editing_data/20210729
mv combined_editing_data.txt $my_dir/pipeline_output/editing_data/20210729
mv *summary.txt $my_dir/pipeline_output/editing_data/20210729
mkdir $my_dir/pipeline_output/variants/20210729
mv *.txt $my_dir/pipeline_output/variants/20210729
mkdir $my_dir/pipeline_output/cigar_counts/20210729
mv ../../cigar_counts/*.txt $my_dir/pipeline_output/cigar_counts/20210729
cd ..
rm *.txt
cd ..
rmdir cigar_counts
rmdir variant_counts/final
rm variant_counts/*.txt
rmdir variant_counts
cd ../..


##### 20210915

cigar_comparisons="VHLx3a_HDRL+VHLx3a_neg,VHLx3a_HDRL+VHLx3ar1_D6,VHLx3ar1_D6+VHLx3ar1_D13,VHLx3ar1_D6+VHLx3ar1_D20,VHLx3ar1_D13+VHLx3ar1_D20,VHLx3a_HDRL+VHLx3a_neg,VHLx3a_HDRL+VHLx3ar2_D6,VHLx3ar2_D6+VHLx3ar2_D13,VHLx3ar2_D6+VHLx3ar2_D20,VHLx3ar2_D13+VHLx3ar2_D20,VHLx3b_HDRL+VHLx3b_neg,VHLx3b_HDRL+VHLx3br1_D6,VHLx3br1_D6+VHLx3br1_D13,VHLx3br1_D6+VHLx3br1_D20,VHLx3br1_D13+VHLx3br1_D20,VHLx3b_HDRL+VHLx3b_neg,VHLx3b_HDRL+VHLx3br2_D6,VHLx3br2_D6+VHLx3br2_D13,VHLx3br2_D6+VHLx3br2_D20,VHLx3br2_D13+VHLx3br2_D20"
amplicon_list="VHLx3a,VHLx3b" #splits on commas
#exp name + 8 sample entries per grouping to be uniform across all runs (some experiments have 5 samples, some 7) must be specified, order must be:  ( experiment + pre + post + HDRLlib + neg + rna + post2 + post3 + rna2 , ... )
exp_groupings="VHLx3arL41+VHLx3ar1_D6+VHLx3ar1_D20+VHLx3a_HDRL+VHLx3a_neg+VHLx3ar1_D13+VHLx3ar1_D13+VHLx3ar1_D13+VHLx3ar1_D13,VHLx3arL42+VHLx3ar2_D6+VHLx3ar2_D20+VHLx3a_HDRL+VHLx3a_neg+VHLx3ar2_D13+VHLx3ar2_D13+VHLx3ar2_D13+VHLx3ar2_D13,VHLx3brL41+VHLx3br1_D6+VHLx3br1_D20+VHLx3b_HDRL+VHLx3b_neg+VHLx3br1_D13+VHLx3br1_D13+VHLx3br1_D13+VHLx3br1_D13,VHLx3brL42+VHLx3br2_D6+VHLx3br2_D20+VHLx3b_HDRL+VHLx3b_neg+VHLx3br2_D13+VHLx3br2_D13+VHLx3br2_D13+VHLx3br2_D13"

mkdir $my_dir/pipeline_output/sam/20210915
cd $my_dir/fastq/20210915

mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns
module load GCC/5.4.0-2.26
module load EMBOSS

mkdir sam
python ~/home/users/findlag/bin/run_needle_to_sam_sbatch_pipeline.py $ref_folder
sbatch run_needle_to_sam_sbatch.sh 

mv sam/* $my_dir/pipeline_output/sam/20210915
cd $my_dir/pipeline_output/sam/20210915

mkdir cigar_counts
python ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -20 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt

mkdir variant_counts
python /camp/lab/findlayg/home/users/findlag/bin/20220826_sam_to_edits_v9s_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold
cd variant_counts
mkdir final
python ~/home/users/findlag/bin/220407_annotate_variants_pipeline.py $amplicon_list $ref_folder $editing_info_file $orientation $CCDS $amplicon_coords $cadd_dir $cadd_version $clinvar_file 
cd final
head -10 *summary.txt >> combined_editing_data.txt

mkdir $my_dir/pipeline_output/editing_data/20210915
mv combined_editing_data.txt $my_dir/pipeline_output/editing_data/20210915
mv *summary.txt $my_dir/pipeline_output/editing_data/20210915
mkdir $my_dir/pipeline_output/variants/20210915
mv *.txt $my_dir/pipeline_output/variants/20210915
mkdir $my_dir/pipeline_output/cigar_counts/20210915
mv ../../cigar_counts/*.txt $my_dir/pipeline_output/cigar_counts/20210915
cd ../..
rmdir cigar_counts
rmdir variant_counts/final
rm variant_counts/*.txt
rmdir variant_counts


##### 20220511

cigar_comparisons="VHLx3b_HDRL+VHLx3b_neg,VHLx3b_HDRL+VHLx3brL44_D6,VHLx3brL44_D6+VHLx3brL44_D13,VHLx3brL44_D6+VHLx3brL44_D20,VHLx3brL44_D13+VHLx3brL44_D20,VHLx3b_HDRL+VHLx3brL43_D6,VHLx3brL43_D6+VHLx3brL43_D13,VHLx3brL43_D6+VHLx3brL43_D20,VHLx3brL43_D13+VHLx3brL43_D20,VHLx2_HDRL+VHLx2_neg,VHLx2_HDRL+VHLx2rL43_D6,VHLx2rL43_D6+VHLx2rL43_D13,VHLx2rL43_D6+VHLx2rL43_D20,VHLx2rL43_D13+VHLx2rL43_D20,VHLx2_HDRL+VHLx2_neg,VHLx2_HDRL+VHLx2rL44_D6,VHLx2rL44_D6+VHLx2rL44_D13,VHLx2rL44_D6+VHLx2rL44_D20,VHLx2rL44_D13+VHLx2rL44_D20,VHLx3a_HDRL+VHLx3a_neg,VHLx3a_HDRL+VHLx3arL43_D6,VHLx3arL43_D6+VHLx3arL43_D13,VHLx3arL43_D6+VHLx3arL43_D20,VHLx3arL43_D13+VHLx3arL43_D20,VHLx3a_HDRL+VHLx3a_neg,VHLx3a_HDRL+VHLx3arL44_D6,VHLx3arL44_D6+VHLx3arL44_D13,VHLx3arL44_D6+VHLx3arL44_D20,VHLx3arL44_D13+VHLx3arL44_D20"
amplicon_list="VHLx2,VHLx3a,VHLx3b" #splits on commas
#exp name + 8 sample entries per grouping to be uniform across all runs (some experiments have 5 samples, some 7) must be specified, order must be:  ( experiment + pre + post + HDRLlib + neg + rna + post2 + post3 + rna2 , ... )
exp_groupings="VHLx2rL43+VHLx2rL43_D6+VHLx2rL43_D20+VHLx2_HDRL+VHLx2_neg+VHLx2rL43_RNAd6gDNA+VHLx2rL43_D13+VHLx2rL43_D13+VHLx2rL43_RNAd20gDNA,VHLx2rL44+VHLx2rL44_D6+VHLx2rL44_D20+VHLx2_HDRL+VHLx2_neg+VHLx2rL44_RNAd6gDNA+VHLx2rL44_D13+VHLx2rL44_D13+VHLx2rL44_RNAd20gDNA,VHLx3arL43+VHLx3arL43_D6+VHLx3arL43_D20+VHLx3a_HDRL+VHLx3a_neg+VHLx3arL43_RNAd6gDNA+VHLx3arL43_D13+VHLx3arL43_D13+VHLx3arL43_RNAd20gDNA,VHLx3arL44+VHLx3arL44_D6+VHLx3arL44_D20+VHLx3a_HDRL+VHLx3a_neg+VHLx3arL44_RNAd6gDNA+VHLx3arL44_D13+VHLx3arL44_D13+VHLx3arL44_RNAd20gDNA,VHLx3brL43+VHLx3brL43_D6+VHLx3brL43_D20+VHLx3b_HDRL+VHLx3b_neg+VHLx3brL43_RNAd6gDNA+VHLx3brL43_D13+VHLx3brL43_D13+VHLx3brL43_RNAd20gDNA,VHLx3brL44+VHLx3brL44_D6+VHLx3brL44_D20+VHLx3b_HDRL+VHLx3b_neg+VHLx3brL44_RNAd6gDNA+VHLx3brL44_D13+VHLx3brL44_D13+VHLx3brL44_RNAd20gDNA"

mkdir $my_dir/pipeline_output/sam/20220511
cd $my_dir/fastq/20220511

mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh

cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns

python ~/home/users/findlag/bin/run_cDNA_to_gDNA_SGE_pipeline.py $cDNA_info_file
sh run_cDNA_to_gDNA.sh

module load GCC/5.4.0-2.26
module load EMBOSS

mkdir sam
python ~/home/users/findlag/bin/run_needle_to_sam_sbatch_pipeline.py $ref_folder
sbatch run_needle_to_sam_sbatch.sh 

mv sam/* $my_dir/pipeline_output/sam/20220511
cd $my_dir/pipeline_output/sam/20220511

mkdir cigar_counts
python ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -20 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt

mkdir variant_counts
python /camp/lab/findlayg/home/users/findlag/bin/20220826_sam_to_edits_v9s_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold
cd variant_counts
mkdir final
python ~/home/users/findlag/bin/220407_annotate_variants_pipeline.py $amplicon_list $ref_folder $editing_info_file $orientation $CCDS $amplicon_coords $cadd_dir $cadd_version $clinvar_file 
cd final
head -10 *summary.txt >> combined_editing_data.txt

mkdir $my_dir/pipeline_output/editing_data/20220511
mv combined_editing_data.txt $my_dir/pipeline_output/editing_data/20220511
mv *summary.txt $my_dir/pipeline_output/editing_data/20220511
mkdir $my_dir/pipeline_output/variants/20220511
mv *.txt $my_dir/pipeline_output/variants/20220511
mkdir $my_dir/pipeline_output/cigar_counts/20220511
mv ../../cigar_counts/*.txt $my_dir/pipeline_output/cigar_counts/20220511

cd ../..
rmdir cigar_counts
rmdir variant_counts/final
rm variant_counts/*.txt
rmdir variant_counts


#20220822

cigar_comparisons="VHLx1a2_HDRL+VHLx1a2_neg,VHLx1a2_HDRL+VHLx1a2rL41_D6,VHLx1a2rL41_D6+VHLx1a2rL41_D13,VHLx1a2rL41_D6+VHLx1a2rL41_D20,VHLx1a2rL41_D13+VHLx1a2rL41_D20,VHLx1a2_HDRL+VHLx1a2rL42_D6,VHLx1a2rL42_D6+VHLx1a2rL42_D13,VHLx1a2rL42_D6+VHLx1a2rL42_D20,VHLx1a2rL42_D13+VHLx1a2rL42_D20,VHLx1b_HDRL+VHLx1b_neg,VHLx1b_HDRL+VHLx1brL43_D6,VHLx1brL43_D6+VHLx1brL43_D13,VHLx1brL43_D6+VHLx1brL43_D20,VHLx1brL43_D13+VHLx1brL43_D20,VHLx1b_HDRL+VHLx1brL44_D6,VHLx1brL44_D6+VHLx1brL44_D13,VHLx1brL44_D6+VHLx1brL44_D20,VHLx1brL44_D13+VHLx1brL44_D20,VHLx1c2_HDRL+VHLx1c2_neg,VHLx1c2_HDRL+VHLx1c2rL41_D6,VHLx1c2rL41_D6+VHLx1c2rL41_D13,VHLx1c2rL41_D6+VHLx1c2rL41_D20,VHLx1c2rL41_D13+VHLx1c2rL41_D20,VHLx1c2_HDRL+VHLx1c2rL42_D6,VHLx1c2rL42_D6+VHLx1c2rL42_D13,VHLx1c2rL42_D6+VHLx1c2rL42_D20,VHLx1c2rL42_D13+VHLx1c2rL42_D20,VHLx2rL4H_HDRL+VHLx2rL4H_neg,VHLx2rL4H_HDRL+VHLx2rL4H1_D6,VHLx2rL4H1_D6+VHLx2rL4H1_D13,VHLx2rL4H1_D6+VHLx2rL4H1_D20,VHLx2rL4H1_D13+VHLx2rL4H1_D20,VHLx2rL4H_HDRL+VHLx2rL4H2_D6,VHLx2rL4H2_D6+VHLx2rL4H2_D13,VHLx2rL4H2_D6+VHLx2rL4H2_D20,VHLx2rL4H2_D13+VHLx2rL4H2_D20,VHLx3arL4H_HDRL+VHLx3arL4H_neg,VHLx3arL4H_HDRL+VHLx3arL4H1_D6,VHLx3arL4H1_D6+VHLx3arL4H1_D13,VHLx3arL4H1_D6+VHLx3arL4H1_D20,VHLx3arL4H1_D13+VHLx3arL4H1_D20,VHLx3arL4H_HDRL+VHLx3arL4H2_D6,VHLx3arL4H2_D6+VHLx3arL4H2_D13,VHLx3arL4H2_D6+VHLx3arL4H2_D20,VHLx3arL4H2_D13+VHLx3arL4H2_D20"
amplicon_list="VHLx1c2,VHLx1b,VHLx1a2,VHLx2,VHLx3a" #splits on commas
#exp name + 8 sample entries per grouping to be uniform across all runs (some experiments have 5 samples, some 7) must be specified, order must be:  ( experiment + pre + post + HDRLlib + neg + rna + post2 + post3 + rna2 , ... )
exp_groupings="VHLx1a2rL41+VHLx1a2rL41_D6+VHLx1a2rL41_D20+VHLx1a2_HDRL+VHLx1a2_neg+VHLx1a2rL41_D6RNAgDNA+VHLx1a2rL41_D13+VHLx1a2rL41_D13+VHLx1a2rL41_D20RNAgDNA,VHLx1a2rL42+VHLx1a2rL42_D6+VHLx1a2rL42_D20+VHLx1a2_HDRL+VHLx1a2_neg+VHLx1a2rL42_D6RNAgDNA+VHLx1a2rL42_D13+VHLx1a2rL42_D13+VHLx1a2rL42_D20RNAgDNA,VHLx1c2rL41+VHLx1c2rL41_D6+VHLx1c2rL41_D20+VHLx1c2_HDRL+VHLx1c2_neg+VHLx1c2rL41_D13+VHLx1c2rL41_D13+VHLx1c2rL41_D13+VHLx1c2rL41_D13,VHLx1c2rL42+VHLx1c2rL42_D6+VHLx1c2rL42_D20+VHLx1c2_HDRL+VHLx1c2_neg+VHLx1c2rL42_D13+VHLx1c2rL42_D13+VHLx1c2rL42_D13+VHLx1c2rL42_D13,VHLx2rL4H1+VHLx2rL4H1_D6+VHLx2rL4H1_D20+VHLx2rL4H_HDRL+VHLx2rL4H_neg+VHLx2rL4H1_D13+VHLx2rL4H1_D13+VHLx2rL4H1_D13+VHLx2rL4H1_D13,VHLx2rL4H2+VHLx2rL4H2_D6+VHLx2rL4H2_D20+VHLx2rL4H_HDRL+VHLx2rL4H_neg+VHLx2rL4H2_D13+VHLx2rL4H2_D13+VHLx2rL4H2_D13+VHLx2rL4H2_D13,VHLx3arL4H1+VHLx3arL4H1_D6+VHLx3arL4H1_D20+VHLx3arL4H_HDRL+VHLx3arL4H_neg+VHLx3arL4H1_D13+VHLx3arL4H1_D13+VHLx3arL4H1_D13+VHLx3arL4H1_D13,VHLx3arL4H2+VHLx3arL4H2_D6+VHLx3arL4H2_D20+VHLx3arL4H_HDRL+VHLx3arL4H_neg+VHLx3arL4H2_D13+VHLx3arL4H2_D13+VHLx3arL4H2_D13+VHLx3arL4H2_D13,VHLx1brL43+VHLx1brL43_D6+VHLx1brL43_D20+VHLx1b_HDRL+VHLx1b_neg+VHLx1brL43_D6RNAgDNA+VHLx1brL43_D13+VHLx1brL43_D13+VHLx1brL43_D20RNAgDNA,VHLx1brL44+VHLx1brL44_D6+VHLx1brL44_D20+VHLx1b_HDRL+VHLx1b_neg+VHLx1brL44_D6RNAgDNA+VHLx1brL44_D13+VHLx1brL44_D13+VHLx1brL44_D20RNAgDNA"

mkdir $my_dir/pipeline_output/sam/20220822
cd $my_dir/fastq/20220822

mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #uses custom remove_n_bases.py script - this step seems only necessary for NextSeq runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns

#The primer used to amplify the VHLx1c2 library was short 5 bp (compared to primer used on VHLx1c2 gDNA samples). Therefore, this is a script to computationally pad the VHLx1c2 library reads for compatibility with downstream scripts relying on exact sequence matching between samples.
python /camp/lab/findlayg/home/users/bucklem/_Scripts/_VHL/Analysis_pipelines/pad_to_gDNA_VHLx1c-hdrl.py VHLx1c2-HDRL.merged.fastq VHLx1c2-HDRL_temp.merged.fastq GAGAT
mv VHLx1c2-HDRL_temp.merged.fastq VHLx1c2-HDRL.merged.fastq

python ~/home/users/findlag/bin/run_cDNA_to_gDNA_SGE_pipeline.py $cDNA_info_file #manually edit output script - permission issue being in greg's directory
sh run_cDNA_to_gDNA.sh

module load GCC/5.4.0-2.26
module load EMBOSS

mkdir sam
python ~/home/users/findlag/bin/run_needle_to_sam_sbatch_pipeline.py $ref_folder
sbatch run_needle_to_sam_sbatch.sh 

mv sam/* $my_dir/pipeline_output/sam/20220822

cd $my_dir/pipeline_output/sam/20220822

mkdir cigar_counts
python ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -20 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt

mkdir variant_counts
python /camp/lab/findlayg/home/users/findlag/bin/20220826_sam_to_edits_v9s_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold
cd variant_counts
mkdir final
python ~/home/users/findlag/bin/220407_annotate_variants_pipeline.py $amplicon_list $ref_folder $editing_info_file $orientation $CCDS $amplicon_coords $cadd_dir $cadd_version $clinvar_file 
cd final
head -10 *summary.txt >> combined_editing_data.txt

mkdir $my_dir/pipeline_output/editing_data/20220822
mv combined_editing_data.txt $my_dir/pipeline_output/editing_data/20220822
mv *summary.txt $my_dir/pipeline_output/editing_data/20220822
mkdir $my_dir/pipeline_output/variants/20220822
mv *.txt $my_dir/pipeline_output/variants/20220822
mkdir $my_dir/pipeline_output/cigar_counts/20220822
mv ../../cigar_counts/*.txt $my_dir/pipeline_output/cigar_counts/20220822
cd ../..
rmdir cigar_counts
rmdir variant_counts/final
rm variant_counts/*.txt
rmdir variant_counts

#20221102


cigar_comparisons="VHLx2rLD1_D6+VHLx2rLD1_D13,VHLx2rLD2_D6+VHLx2rLD2_D13,VHLx2rLDB1_D6+VHLx2rLDB1_D13,VHLx2rLDB2_D6+VHLx2rLDB2_D13,VHLx2rLDH1_D6+VHLx2rLDH1_D13,VHLx2rLDH2_D6+VHLx2rLDH2_D13,VHLx1p_neg+VHLx1p_HDRL,VHLx1p_HDRL+VHLx1prL41_D6,VHLx1p_HDRL+VHLx1prL42_D6,VHLx1prL41_D6+VHLx1prL41_D13,VHLx1prL41_D6+VHLx1prL41_D20,VHLx1prL41_D13+VHLx1prL41_D20,VHLx1prL42_D6+VHLx1prL42_D13,VHLx1prL42_D6+VHLx1prL42_D20,VHLx1prL42_D13+VHLx1prL42_D20"
amplicon_list="VHLx1p,VHLx2" #splits on commas
#exp name + 8 sample entries per grouping to be uniform across all runs (some experiments have 5 samples, some 7) must be specified, order must be:  ( experiment + pre + post + HDRL + neg + rna + post2 + post3 + rna2 )
exp_groupings="VHLx1prL41+VHLx1prL41_D6+VHLx1prL41_D20+VHLx1p_HDRL+VHLx1p_neg+VHLx1prL41_D13+VHLx1prL41_D13+VHLx1prL41_D13+VHLx1prL41_D13,VHLx1prL42+VHLx1prL42_D6+VHLx1prL42_D20+VHLx1p_HDRL+VHLx1p_neg+VHLx1prL42_D13+VHLx1prL42_D13+VHLx1prL42_D13+VHLx1prL42_D13"

mkdir $my_dir/pipeline_output/sam/20221102
cd $my_dir/fastq/20221102

mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns

module load GCC/5.4.0-2.26
module load EMBOSS

mkdir sam
python ~/home/users/findlag/bin/run_needle_to_sam_sbatch_pipeline.py $ref_folder

sbatch run_needle_to_sam_sbatch.sh 

mv sam/* $my_dir/pipeline_output/sam/20221102

cd $my_dir/pipeline_output/sam/20221102
mkdir cigar_counts
python ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -20 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt

mkdir variant_counts
python /camp/lab/findlayg/home/users/findlag/bin/20220826_sam_to_edits_v9s_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold
cd variant_counts
mkdir final
python ~/home/users/findlag/bin/220407_annotate_variants_pipeline.py $amplicon_list $ref_folder $editing_info_file $orientation $CCDS $amplicon_coords $cadd_dir $cadd_version $clinvar_file 
cd final
head -10 *summary.txt >> combined_editing_data.txt

mkdir $my_dir/pipeline_output/editing_data/20221102
mv combined_editing_data.txt $my_dir/pipeline_output/editing_data/20221102
mv *summary.txt $my_dir/pipeline_output/editing_data/20221102
mkdir $my_dir/pipeline_output/variants/20221102
mv *.txt $my_dir/pipeline_output/variants/20221102
mkdir $my_dir/pipeline_output/cigar_counts/20221102
mv ../../cigar_counts/*.txt $my_dir/pipeline_output/cigar_counts/20221102
cd ../..
rmdir cigar_counts
rmdir variant_counts/final
rm variant_counts/*.txt
rmdir variant_counts

#####END, subsequent analyses performed in R using output files in the 'variants' folder: 

