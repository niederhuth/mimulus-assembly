#! /usr/bin/perl

# rename_fastas_with_new_identifier_v2.pl

# 23 March 2017

# Kevin Childs

# Updated 21 June 2021

# Given a fasta file with awkward identifiers, make a new fasta file with
# new identifiers.  A file with the new-to-old identifiers is also required.



use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use strict;
use warnings;

my $usage = "\n$0\n    --input_fasta  input_fastq_with_poor_ids\n" .
                  "    --output_fasta   output_fastq_with_better_ids\n" .
                  "    --new_ids_list   file_with_new_and_old_gene_ids\n" .
                  "    [--help]\n\n";

my ($input_fasta, $output_fasta, $new_ids_list);
my $help;

Getopt::Long::GetOptions( "input_fasta=s" => \$input_fasta,
                          "output_fasta=s" => \$output_fasta,
                          "new_ids_list=s" => \$new_ids_list,
                          "help" => \$help) || die;

if (defined($help)) {
    print $usage;
    exit;
}

if (!defined($input_fasta) ||  !(-e $input_fasta) ||
    !defined($output_fasta) ||  (-e $output_fasta) ||
    !defined($new_ids_list) ||  !(-e $new_ids_list)) {
    die $usage;
}

my %old_to_new_ids;

open IN, $new_ids_list or die "\nUnable to open $new_ids_list for reading.\n\n";

while (my $line = <IN>) {

    chomp $line;
    my ($new_id, $old_id) = split "\t", $line;

    print "new $new_id\told $old_id\n";

    $old_to_new_ids{$old_id} = $new_id;

}


my $input = Bio::SeqIO->new(-file => $input_fasta, -format => "fasta");
my $output = Bio::SeqIO->new(-file => ">$output_fasta", -format => "fasta");
while (my $seq_obj = $input->next_seq()) {
    my $id = $seq_obj->id();

    if (exists($old_to_new_ids{$id})) {
	$seq_obj->id($old_to_new_ids{$id});
	$output->write_seq($seq_obj);
    }
    else {
	die "\nThe identifier, $id, was not found in the decoder list.\n\n";
    }
}

exit;

