#I really hate using web interfaces for my bioinformatics. Especially when I do all
#my real work on a remote server. Also, I don't really have the desire to maintain an
#up to date full BLAST database locally for random projects. This script is designed to
#let me do remote Blast searches at NCBI and convert and filter the results for easy input 
#into downstream applications, such as building gene trees. I am sure this is also not 
#the prettiest code, but it works. ~Chad Niederhuth

#Requires argparse and biopython
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

#Function for running NCBIWWW.qblast
def my_blast(blast_type, database, input, output, results=100, evalue=0.0001, identity=None):
	#Read input fasta with SeqIO
	print("Reading input fasta")
	query = SeqIO.read(input, format="fasta")
	#Run BLAST
	print("Running NCBI BLAST")
	print("This may take a while")
	blast_result = NCBIWWW.qblast(blast_type, database, query.format("fasta"),
		hitlist_size=results, expect=evalue, perc_ident=identity)
	print("BLAST complete")
	#Write blast_results to xml output
	print("Converting results to XML")
	with open(output, "w") as out_xml:
		out_xml.write(blast_result.read())
	#Close files
	blast_result.close()
	out_xml.close()

#Function for downloading genbank records
def download_gb(db, input, output):
	#Read the XML file and create an iterator of each record
	in_xml = open(input, "r")
	blast_record = NCBIXML.read(in_xml)
	in_xml.close()
	#Download genbank records and write to single genbank output file
	#Why genbank instead of fasta? Because genbank has additional information that allows us to 
	#filter the results. For instance if you want to limit it to specific taxonomic groups
	print("Downloading genbank results")
	#Open output file
	id_list = []
	with open(output, "w") as out_gb:
		#Iterate over blast records 
		for alignment in blast_record.alignments:
			id_list += [str(alignment.accession)]
			#For each blast record, download genbank record for that accession
		gb_result = Entrez.efetch(db=db, id=id_list, rettype="gb", retmode="text")
		out_gb.write(gb_result.read())
	#Close the output genbank file
	out_gb.close()

#Function for filtering and converting genbank to fasta
def gb2fa(input, output, include=[], exclude=[]):
	#Filter genbank and convert to fasta
	print("Converting to fasta")
	with open(input, "r") as input_gb, open(output, "w") as out_fa:
		#Set count to 0
		count=0
		#Read in the genbank records
		sequences = SeqIO.parse(input_gb, "genbank")
		#Iterate over each genbank record
		for record in sequences:
			#Exclude records that are part of the --exclude listed taxonomic groups
			if not any(b in record.annotations.get('taxonomy') for b in exclude):
				#Check if inlcude option is not empty
				if include:
					#Include only those records that are part of the --include listed taxonomic groups
					if any(a in record.annotations.get('taxonomy') for a in include):
						#Use SeqIO to convert genbank to fasta
						count = count + SeqIO.write(record, out_fa, "fasta")
				else:
					count = count + SeqIO.write(record, out_fa, "fasta")
			else:
				#If a record is not part of --include list or is part of --exclude list, 
				#report that you are skipping it. These are tab-delimited so to make it easier to extract
				#in case you want to check which records are skipped.
				print("Skipping" + '\t' + record.id)
	#Close files
	input_gb.close()
	out_fa.close()
	#Report how many records were converted to fasta
	print("Converted %i records" % count)

if __name__ == "__main__":
	#Set up arguments options and parsing
	#My first time trying this, so I'm sure its ugly
	parser = argparse.ArgumentParser(description='This script is for running Blast \
		remotely on NCBI, downloading, filtering, and converting the results. There \
		are three key steps, each of which can be run independently. Step 1: Run \
		Blast with your input sequence remotely on NCBI and download the results in \
		XML format. Step 2: Download genbank records for each Blast result. Step 3: \
		Filter genbank records and convert to fasta. Right now filtering in the 3rd \
		step can only be done based on the taxonomic lineage of that sequence. There \
		are two options for filtering. You can include any sequence within a lineage \
		(e.g. species, genus, etc) or you can exclude any sequence in a lineage. This \
		does require you to know how those lineages you want to include or exclude are \
		spelled (including capitalization) in genbank. However this can be useful if say \
		you want to only include results for a specific species (e.g. Vitis vinifera) or \
		a lineage like green plants (Viridiplantae). More options can easily be added by \
		modifying this script.')

	#Add a required argument for the output data file prefix
	parser.add_argument('--output-prefix', required=True, dest='output_prefix', 
		help='Prefix for output files')

	#Add an optional argument for the input data file name and open in 'read' mode
	#This is required to run Blast!
	parser.add_argument('--input', type=argparse.FileType('r'), 
		help='Input fasta file of sequence to Blast. This is only required if \
		running Blast. As of now, only run this with a single sequence in the file. \
		I have not tested it with multiple sequences and latter steps are not set up \
		to properly parse the results.')
	#Add an optional argument for blast progam
	#i.e. blastp, blastn, etc
	parser.add_argument('--bp', '--blast-program', dest='blast_program', 
		help='Blast program used. This is required to run Blast and download genbank files')
	#Add an optional argument for blast database
	parser.add_argument('--db', 
		help='Which Blast database to use, e.g. "nr" for blastp')
	#Add an optional argument for number of results to return, default = 100
	parser.add_argument('--rn', '--result-number', type=int, dest='result_number', default=100,
		help='Number of Blast results to return. The default is 100, but I would set this \
		much higher for most applications.')
	#Add an optional argument for number for evalue cutoff
	parser.add_argument('--evalue',  default=0.0001,
		help='E-value cutoff. The defualt is 0.0001')
	#Add an optional argument for number for evalue cutoff
	parser.add_argument('--identity',  type=int, default=None,
		help='Percent identity. The default is None')
	#Add an optional argument for included species/lineages
	#Use nargs="*" so that you can list multiple arguments
	parser.add_argument('--include', nargs='*',
		help='Taxonomic groups to include in the fasta file during conversion. \
		This can be a space delimited list, e.g "--include Vitis Prunus"')
	#Add an optional argument for excluded species/lineages
	#Use nargs="*" so that you can list multiple arguments
	parser.add_argument('--exclude', nargs='*',
		help='Taxonomic groups to exclude from the fasta file during conversion. \
		This can be a space delimited list, e.g "--exclude Vitis Prunus"')
	#Add your email, used with Entrez efetch
	parser.add_argument('--email',
		help='If you download genbank files, Entrez wants your email. If not provided \
		it will usually still work, but will give an annoying warning message')
	#Add an optional argument to stop after Blast result
	#Will not proceed to genbank and fasta conversion 
	parser.add_argument('--blast-only', action="store_true", dest='blast_only',
		help='If this argument is provided, then only the NCBI Blast will be run')
	#Add an optional argument to skip Blast and download & convert results
	#This assumes you have a blast result in xml format, blast program is still
	#needed, as it informs the Entrez database to download from.
	parser.add_argument('--no-blast', action="store_true", dest='no_blast',
		help='If this argument is provided, then Blast will be skipped and only the \
		genbank download and fasta conversion processes run. This still requires you \
		to provide the Blast program used, as this is used to determine the appropriate \
		Entrez database to use for downloading genbank records. This also assumes you \
		already have Blast results in XML format. These are indicated by providing the \
		required --output-prefix option, using the same prefix that your Blast XML file \
		is named. So if your Blast results are "my_blast.xml", then run with \
		"--output-prefix my_blast".')
	#Add an optional argument to only download genbank results
	#This assumes you have a blast result in xml format, the prefix of which must be supplied.
	#blast program is still needed, as it informs the Entrez database to download from.
	parser.add_argument('--gb-download', action="store_true", dest='gb_download',
		help='If this argument is provided, only the genbank records will be downloaded. \
		This still requires you to provide the Blast program used, as this is used to \
		determine the appropriate Entrez database to use for downloading genbank records. \
		This also assumes you already have Blast results in XML format. These are indicated \
		by providing the required --output-prefix option, using the same prefix that your \
		Blast XML file is named. So if your Blast results are "my_blast.xml", then run with \
		"--output-prefix my_blast".')
	#Add an optional argument to only convert and filter genbank results
	#This assumes you have a genbank records, the prefix of which must be supplied
	parser.add_argument('--convert', action="store_true",
		help='If this argument is provided, only the fasta conversion and filtering will \
		be run. It is not necessary to provide the blast program. It does assume you have \
		genbank records. These are indicated by providing the required --output-prefix option, \
		using the same prefix that your genbank file is named. So if your genbank files are \
		"my_blast.gb", then run with "--output-prefix my_blast".')

	#Parse out the arguments
	args = parser.parse_args()

	#Set blast database. Default is "nr" for protein, nr/nt for nucleotide
	#--db overrides whatever you use as the blast program
	if not args.convert:
		if args.db:
			blast_db=args.db
		#if --db not provided, check blast_program and set to appropriate database
		elif not args.db and args.blast_program in ['blastp','blastx']:
			blast_db='nr'
		elif not args.db and args.blast_program in ['blastn','tblastx','tblastn']:
			blast_db='nr/nt' #I need to double check this!

	#What db to use for Entrez, e.g. protein or nucleotide
	if not args.convert:
		if args.blast_program in ['blastp','blastx']:
			Entrez_db='protein'
		elif args.blast_program in ['blastn','tblastx','tblastn']:
			Entrez_db='nucleotide'

	#Set output files
	file1=args.output_prefix + '.xml'
	file2=args.output_prefix + '.gb'
	file3=args.output_prefix + '.fa'

	#Run Blast if the following arguments are not set
	if not args.no_blast and not args.gb_download and not args.convert:
		my_blast(args.blast_program, blast_db, args.input, file1, results=args.result_number,
			evalue=args.evalue, identity=args.identity)

	#Run download_gb if following arguments are not set
	if not args.blast_only and not args.convert:
		#Set email if provided, will give error message if not provided
		if args.email:
			Entrez.email=args.email
		download_gb(Entrez_db, file1, file2)

	#Run gb2fa if the following arguments are not set
	if not args.blast_only and not args.gb_download:
		#Set lists for inclusion or exclusion based on arguments provided
		if args.include:
			include_list=args.include
		else:
			include_list=[]
		if args.exclude:
			exclude_list=args.exclude
		else:
			exclude_list=[]
		#Filter and Convert
		gb2fa(file2, file3, include=include_list, exclude=exclude_list)

	print("Done")
