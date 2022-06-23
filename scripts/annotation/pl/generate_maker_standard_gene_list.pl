#! /usr/bin/perl

#Original by Kevin Childs 9/20/2013 generate_maker_standard_gene_list.pl
#Modified by Chad Niederhuth 6/1/2022

use Getopt::Long;
use strict;
use warnings;

# We ain't waiting on no stinking print buffer.
$| = 1;

my $usage = "\n$0\n    --input_gff   <path_to_input_gff_file>\n" .
            "    --output_file   <MAKER-standard gene list output file>\n" .
            "    [--help]\n\n";

my ($input_gff, $pfam_results, $pfam_cutoff, $output_file);
my $help;

Getopt::Long::GetOptions( "input_gff=s" => \$input_gff,
                          "output_file=s" => \$output_file,
                          "help" => \$help) || die;

if (defined($help)) {
    print $usage;
    exit;
}
if (!defined($input_gff) ||  !(-e $input_gff) || 
    !defined($output_file) ||  (-e $output_file)) {
    die $usage;
}


open GFF, $input_gff or die "\nUnable to open $input_gff for reading.\n\n";
open OUT, ">$output_file" or die "\nUnable to open $output_file for writing.\n\n";

my %maker_standard_genes;
while (my $line = <GFF>) {
    chomp $line;

    if ($line =~ /#FASTA/) {
	last;
    }
    if ($line =~ /^#/) {
	next;
    }

    my @elems = split "\t", $line;

    if ($elems[2] eq 'gene') {
	if ($elems[8] =~ /ID=([^;]+)/) {
            my $id = $1;
	    print OUT "$id\n";
	}
    }

}
close GFF;
close OUT;

exit;
