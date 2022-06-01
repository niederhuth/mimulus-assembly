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