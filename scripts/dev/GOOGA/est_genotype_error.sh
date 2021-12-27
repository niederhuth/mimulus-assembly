#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=50GB
#SBATCH --job-name est_genotype_error
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
max_error=0.2

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="radseq"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)

#Create bad_marks.txt
cd g_files
touch bad_marks.txt

#Loop over each genotype file and calculate error rates
for i in *_genotype.txt
do
	python ${path1}/dev/GOOGA/hmm_scaff_likelihood.py ${i/_genotype.txt/} >> error_rates_1.txt
done

#Filter samples based on max error rate
awk -v a=${max_error} '{ if ($2<a && $3<a && $4<a) print $1}' error_rates_1.txt >> low_error_samples.txt

#1st linkage map
python ${path1}/dev/GOOGA/hmm_intrascaff_rates.py error_rates_1.txt $(head -1 low_error_samples.txt)_genotype.txt low_error_samples.txt MLE_iscaff_1.txt bad_marks.txt

#Filter markers

#2nd linkage map
python ${path1}/dev/GOOGA/hmm_intrascaff_rates.py error_rates_1.txt $(head -1 low_error_samples.txt)_genotype.txt low_error_samples.txt MLE_iscaff_2.txt bad_marks.txt​

#Loop over each genotype file and recalculate error rates
for i in *_genotype.txt
do
	python ${path1}/dev/GOOGA/hmm_scaff_likelihood.py ${i/_genotype.txt/} >> error_rates_2.txt
done

#Filter samples based on max error rate
awk -v a=${max_error} '{ if ($3<a && $4<a && $5<a) print $1}' error_rates_1.txt >> low_error_samples_2.txt

#3rd linkage map
python ${path1}/dev/GOOGA/hmm_intrascaff_rates.py error_rates_2.txt $(head -1 low_error_samples_2.txt)_genotype.txt low_error_samples_2.txt MLE_iscaff_3.txt bad_marks.txt​

#Something

#genotype each scaffold?
for i in 
do
	python {path1}/dev/GOOGA/genotyp_pp.py ${i} low_error_samples_2.txt all_samples.txt MLE_iscaff_3.txt
done

