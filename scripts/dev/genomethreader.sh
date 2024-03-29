#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name protein2genome
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta=$(ls -l repeatmasker/*.fa.masked | sed s/.*\ //)
proteins=$(ls -l proteins/*fa | sed s/.*\ //)
splice_site_model="undefined"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path2="exonerate"

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Loop over protein sources and run exonerate
for i in ${proteins}
do
	echo "Working on ${i} proteins"
	outdir=$(echo ${i} | sed s/.*\\/// | sed s/\.fa//)
	if [ -d ${outdir} ]
	then
		echo "${outdir} directory present, checking files"
	else
		mkdir ${outdir}
	fi
	if [ -z  ]
	then
		echo "Aligning on ${fasta}"
		gth \
			-genomic ${fasta} \
			-protein ${i} \
            -gff3out \
            -o genomethreader
            #-species ${splice_site_model} \
            #-xmlout            show output in XML format
                   #default: no
            #-md5ids            show MD5 fingerprints as sequence IDs
                   #default: no
            #-skipalignmentout  skip output of spliced alignments
                  # default: no
            #-mincutoffs        show full spliced alignments
            #       i.e., cutoffs mode for leading and terminal bases is MINIMAL
            #       default: no
            #-showintronmaxlen  set the maximum length of a fully shown intron
            #       If set to 0, all introns are shown completely
             #      default: 120
            #-minorflen         set the minimum length of an ORF to be shown
            #       default: 64
            #-startcodon        require than an ORF must begin with a start codon
            #       default: no
            #-finalstopcodon    require that the final ORF must end with a stop codon
             #     default: no
            #-showseqnums       show sequence numbers in output
            #       default: no
            #-pglgentemplate    show genomic template in PGL lines
           #       (switch off for backward compatibility)
           #        default: yes
           # -gs2out            output in old GeneSeqer2 format
            #       default: no
            #-maskpolyatails    mask poly(A) tails in cDNA/EST files
            #       default: no
            #-proteinsmap       specify smap file used for protein files
            #       default: protein
            #-noautoindex       do not create indices automatically
            #       except for the .dna.* files used for the DP.
            #       existence is not tested before an index is actually used!
            #       default: no
            #-createindicesonly stop program flow after the indices have been created
             #      default: no
            #-skipindexcheck    skip index check (in preprocessing phase)
            #       default: no
            #-minmatchlen       specify minimum match length (cDNA matching)
            #       default: 20
            #-seedlength        specify the seed length (cDNA matching)
            #      default: 18
            #-exdrop            specify the Xdrop value for edit distance extension (cDNA
            #       matching)
            #       default: 2
            #-prminmatchlen     specify minimum match length (protein matches)
            #       default: 24
            #-prseedlength      specify seed length (protein matching)
            #       default: 10
            #-prhdist           specify Hamming distance (protein matching)
            #       default: 4
            #-online            run the similarity filter online without using the complete
             #      index (increases runtime)
            #       default: no
            #-inverse           invert query and index in vmatch call
            #       default: no
            #-exact             use exact matches in the similarity filter
             #      default: no
            #-gcmaxgapwidth     set the maximum gap width for global chains
             #      defines approximately the maximum intron length
             #      set to 0 to allow for unlimited length
             #      in order to avoid false-positive exons (lonely exons) at the
             #      sequence ends, it is very important to set this parameter
              #     appropriately!
              #     default: 1000000
            #-gcmincoverage     set the minimum coverage of global chains regarding to the
            #       reference sequence
            #       default: 50
            #-paralogs          compute paralogous genes (different chaining procedure)
            #       default: no
            #-enrichchains      enrich genomic sequence part of global chains with additional
            #       matches
            #       default: no
            #-introncutout      enable the intron cutout technique
            #       default: no
            #-fastdp            use jump table to increase speed of DP calculation
             #      default: no
            #-autointroncutout  set the automatic intron cutout matrix size in megabytes and
            #       enable the automatic intron cutout technique
            #      default: 0
            #-icinitialdelta    set the initial delta used for intron cutouts
            #       default: 50
            #-iciterations      set the number of intron cutout iterations
            #       default: 2
            #-icdeltaincrease   set the delta increase during every iteration
            #       default: 50
            #-icminremintronlen set the minimum remaining intron length for an intron to be
            #       cut out
             #      default: 10
            #-nou12intronmodel  disable the U12-type intron model
            #       default: no
            #-u12donorprob      set the probability for perfect U12-type donor sites
            #       default: 0.99
            #-u12donorprob1mism set the prob. for U12-type donor w. 1 mismatch
            #       default: 0.90
            #-probies           set the initial exon state probability
            #       default: 0.50
            #-probdelgen        set the genomic sequence deletion probability
            #      default: 0.03
            #-identityweight    set the pairs of identical characters weight
             #      default: 2.00
            #-mismatchweight    set the weight for mismatching characters
             #      default: -2.00
            #-undetcharweight   set the weight for undetermined characters
            #       default: 0.00
            #-deletionweight    set the weight for deletions
            #       default: -5.00
            #-dpminexonlen      set the minimum exon length for the DP
            #       default: 5
            #-dpminintronlen    set the minimum intron length for the DP
             #      default: 50
            #-shortexonpenal    set the short exon penalty
             #      default: 100.00
            #-shortintronpenal  set the short intron penalty
            #       default: 100.00
            #-wzerotransition   set the zero transition weights window size
            #       default: 80
           # -wdecreasedoutput  set the decreased output weights window size
             #      default: 80
            #-leadcutoffsmode   set the cutoffs mode for leading bases
             #      can be either RELAXED, STRICT, or MINIMAL
             #      default: RELAXED
            #-termcutoffsmode   set the cutoffs mode for terminal bases
            #       can be either RELAXED, STRICT, or MINIMAL
             #      default: STRICT
            #-cutoffsminexonlen set the cutoffs minimum exon length
            #       default: 5
            #-scoreminexonlen   set the score minimum exon length
            #       default: 50
            #-minaveragessp     set the minimum average splice site prob.
            #       default: 0.50
            #-duplicatecheck    criterion used to check for spliced alignment duplicates,
            #       choose from none|id|desc|seq|both
             #      default: both
            #-minalignmentscore set the minimum alignment score for spliced alignments to be
            #       included into the set of spliced alignments
            #       default: 0.00
            #-maxalignmentscore set the maximum alignment score for spliced alignments to be
            #       included into the set of spliced alignments
            #       default: 1.00
            #-mincoverage       set the minimum coverage for spliced alignments to be
             #      included into the set of spliced alignments
             #      default: 0.00
            #-maxcoverage       set the maximum coverage for spliced alignments to be
            #       included into the set of spliced alignments
            #       default: 9999.99
            #-intermediate      stop after calculation of spliced alignments and output
            #       results in reusable XML format. Do not process this output
             #      yourself, use the ``normal'' XML output instead!
            #       default: no
            #-sortags           sort alternative gene structures according to the weighted
            #      mean of the average exon score and the average splice site
            #       probability
            #      default: no
            #-sortagswf         set the weight factor for the sorting of AGSs
             #      default: 1.00
            #-exondistri        show the exon length distribution
             #      default: no
            #-introndistri      show the intron length distribution
             #      default: no
            #-refseqcovdistri   show the reference sequence coverage distribution
             #      default: no
            #-first             set the maximum number of spliced alignments per genomic DNA
            #       input. Set to 0 for unlimited number.
            #       default: 0
	fi
done

echo "Done"
