#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name=rename_annotations_by_chr
#SBATCH --output=../job_reports/%x-%j.SLURMout

input_gff=
prefix="MgL1"
justify=4

#Get chromosome list
cut -f1 ${input_gff} | grep -v \# | sort | uniq > chr_list

#Loop over chr_list and rename genes for each chromosome
mkdir tmp
cd tmp
cat ../chr_list | while read line
do
	echo "Working on ${line}"
	awk -v a=${line} '$1==a' ../${input_gff} > ${line}.gff
	maker_map_ids \
		--prefix ${prefix}${line}G \
		--justify ${justify} \
		--iterate 1 \
		${line}.gff | awk '{if ($2 ~ /-R/) print $0; else print $0"0"}' | sed s/-R/0./ > ${line}_renamed.gff
done




#

#

maker2eval ${gff}

validate_gtf.pl ${gtf} > maker_std_noTE_with_fasta.gtf_validation

get_general_stats.pl
  -g: Input files are gtf not lists
  -q: Quick load the gtf file.  Do not check them for errors.
  -A: Do not get stats for alternative splices. (Faster)
  -v: Verbose mode
  -h: Display this help message and exit

  # create naming table (there are additional options for naming beyond defaults)
maker_map_ids --prefix BoaCon --justify 5  Bcon_rnd3.all.maker.gff > Bcon_rnd3.all.maker.name.map
# replace names in GFF files
map_gff_ids Bcon_rnd3.all.maker.name.map Bcon_rnd3.all.maker.gff
map_gff_ids Bcon_rnd3.all.maker.name.map Bcon_rnd3.all.maker.noseq.gff
# replace names in FASTA headers
map_fasta_ids Bcon_rnd3.all.maker.name.map Bcon_rnd3.all.maker.transcripts.fasta
map_fasta_ids Bcon_rnd3.all.maker.name.map Bcon_rnd3.all.maker.proteins.fasta


maker_map_ids --prefix MgS1_ --suffix . --iterate 1 --justify 8 S1-v1.gff > test