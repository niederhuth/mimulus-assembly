#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name gene_functions
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/gene-function/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/gene-function/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path3="gene_functions"

#Look for fasta file, there can only be one!
if [ -z ${fasta} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		fasta=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${fasta} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		fasta=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${fasta} found"
	elif ls *.fna >/dev/null 2>&1
	then
		fasta=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${fasta} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${fasta}"
fi

#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#Set various datasets
if [ -z ${transcripts} ]
then
	transcripts=${maker_dir}/${fasta/.fa/}.all.maker.transcripts.fasta
fi
if [ -z ${proteins} ]
then
	proteins=${maker_dir}/${fasta/.fa/}.all.maker.proteins.fasta
fi
if [ -z ${gff} ]
then
	gff=${maker_dir}/${fasta/.fa/}.all.gff
fi

#Run InteProScan
echo ""
interproscan.sh \
	-appl pfam \
	-dp \
	-f TSV \
	-goterms \
	-iprlookup \
	-pa \
	-t p \
	-cpu ${threads} \
	-i proteins_noTE.fasta \
	-o proteins.final.fasta.iprscan


 -b,--output-file-base <OUTPUT-FILE-BASE>                  Optional, base output filename (relative or absolute path).
                                                           Note that this option, the --output-dir (-d) option and the
                                                           --outfile (-o) option are mutually exclusive.  The
                                                           appropriate file extension for the output format(s) will be
                                                           appended automatically. By default the input file path/name
                                                           will be used.
 -cpu,--cpu <CPU>                                          Optional, number of cores for inteproscan.
 -crid,--clusterrunid <CLUSTER-RUN-ID>                     Optional, switch to specify the Project name for this i5 run.
 -d,--output-dir <OUTPUT-DIR>                              Optional, output directory.  Note that this option, the
                                                           --outfile (-o) option and the --output-file-base (-b) option
                                                           are mutually exclusive. The output filename(s) are the same
                                                           as the input filename, with the appropriate file extension(s)
                                                           for the output format(s) appended automatically .
 -dp,--disable-precalc                                     Optional.  Disables use of the precalculated match lookup
                                                           service.  All match calculations will be run locally.
 -dra,--disable-residue-annot                              Optional, excludes sites from the XML, JSON output
 -etra,--enable-tsv-residue-annot                          Optional, includes sites in TSV output
 -exclappl,--excl-applications <EXC-ANALYSES>              Optional, comma separated list of analyses you want to
                                                           exclude.
 -f,--formats <OUTPUT-FORMATS>                             Optional, case-insensitive, comma separated list of output
                                                           formats. Supported formats are TSV, XML, JSON, GFF3, HTML and
                                                           SVG. Default for protein sequences are TSV, XML and GFF3, or
                                                           for nucleotide sequences GFF3 and XML.
 -goterms,--goterms                                        Optional, switch on lookup of corresponding Gene Ontology
                                                           annotation (IMPLIES -iprlookup option)
 -help,--help                                              Optional, display help information
 -hm,--highmem                                             Optional, switch on the creation of a high memory worker.
                                                           Please note normal and high mem workers share the same Spring
                                                           configuration file.
 -i,--input <INPUT-FILE-PATH>                              Optional, path to fasta file that should be loaded on Master
                                                           startup. Alternatively, in CONVERT mode, the InterProScan 5
                                                           XML file to convert.
 -incldepappl,--incl-dep-applications <INC-DEP-ANALYSES>   Optional, comma separated list of deprecated analyses that
                                                           you want included.  If this option is not set, deprecated
                                                           analyses will not run.
 -iprlookup,--iprlookup                                    Also include lookup of corresponding InterPro annotation in
                                                           the TSV and GFF3 output formats.
 -m,--mode <MODE-NAME>                                     Optional, the mode in which InterProScan is being run, the
                                                           default mode is standalone. Must be one of: master, worker,
                                                           distributed_worker, highmem_worker, standalone,
                                                           distributed_master, cluster, singleseq, installer,
                                                           empty_installer, convert, monitor.
 -mastermaxlife,--mastermaxlife <MASTER-MAXLIFE>           The maximum lifetime of the Master.
 -masteruri,--masteruri <MASTER-URI>                       The TCP URI of the Master.
 -ms,--minsize <MINIMUM-SIZE>                              Optional, minimum nucleotide size of ORF to report. Will only
                                                           be considered if n is specified as a sequence type. Please be
                                                           aware of the fact that if you specify a too short value it
                                                           might be that the analysis takes a very long time!
 -o,--outfile <EXPLICIT_OUTPUT_FILENAME>                   Optional explicit output file name (relative or absolute
                                                           path).  Note that this option, the --output-dir (-d) option
                                                           and the --output-file-base (-b) option are mutually
                                                           exclusive. If this option is given, you MUST specify a single
                                                           output format using the -f option.  The output file name will
                                                           not be modified. Note that specifying an output file name
                                                           using this option OVERWRITES ANY EXISTING FILE.
 -p,--priority <JMS-PRIORITY>                              Minimum message priority that the worker will accept (0 low
                                                           -> 9 high).
 -pa,--pathways                                            Optional, switch on lookup of corresponding Pathway
                                                           annotation (IMPLIES -iprlookup option)
 -sct,--sequencecount <SEQUENCE_COUNT>                     The total number of input sequences in master.
 -t,--seqtype <SEQUENCE-TYPE>                              Optional, the type of the input sequences (dna/rna (n) or
                                                           protein (p)).  The default sequence type is protein.
 -T,--tempdir <TEMP-DIR>                                   Optional, specify temporary file directory (relative or
                                                           absolute path). The default location is temp/.
 -td,--tempdirname <TEMP-DIR-NAME>                         Optional, used to start up a worker with the correct
                                                           temporary directory.
 -tier1,--tier1 <TIER>                                     Optional, switch to indicate the high memory worker is a
                                                           child of the master.
 -u,--userdir <USER_DIRECTORY>                             The base directory for results (if absolute paths not
                                                           specified)
 -verbose,--verbose                                        Optional, display more verbose log output
 -version,--version                                        Optional, display version number
 -vl,--verbose-level <VERBOSE-LEVEL>                       Optional, display verbose log output at level specified.
 -vtsv,--output-tsv-version                                Optional, includes a TSV version file alon

makeblastdb \
	-in TAIR10 \
	-input_type fasta \
	-dbtype prot \
	-title "TAIR10DB" \
	-parse_seqids \
	-out TAIR10DB
blastp \
	-query proteins_noTE.fasta \
	-db TAIR10DB \
	-evalue 1e-6 \
	-max_hsps 1 \
	-max_target_seqs 5 \
	-num_threads 15 \
	-out Sreb_TAIR10_blast.out

perl -e  'while (my $line = <>){ my @elems = split "\t", $line; if($elems[2] ne "")\
{print "$elems[0]\t$elems[2]\n"}}' TAIR10_functional_descriptions > TAIR10_short_functional_descriptions.txt

protein=proteins_noTE.fasta
annot=TAIR10_short_functional_descriptions.txt
blast=blast.out
pfam=prot_domains.out

create_functional_annotation_file_2020.pl \
	--protein_fasta $protein \
	--model_annot $annot \
	--model_blast $blast \
	--pfam_results_file $pfam \ 
	--max_hits 5 \
	--output prot_func_desc_list.txt