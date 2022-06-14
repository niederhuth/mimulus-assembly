#! /usr/bin/perl

# rename_gff_with_new_identifier_v2.pl

# 23 March 2017

# Kevin Childs

# Given a maker gff3 file with awkward maker identifiers, make a new gff3 file with
# new identifiers.  A file with the new-to-old identifiers is also required.

# This new version will give transcripts and proteins an extension that will correspond to their isoform number.
# This has not been an issue with us to this point, but we may retain multiple isoforms in the future, and this
# change will allow us to keep track of those.


use Getopt::Long;
use strict;
use warnings;

my $usage = "\n$0\n    --input_gff  input_gff_with_maker_ids\n" .
                  "    --output_gff   output_gff_with_new_ids\n" .
                  "    --new_ids_list   file_with_new_and_old_gene_ids\n" .
                  "    [--help]\n\n";

my ($input_gff, $output_gff, $new_ids_list);
my $help;

Getopt::Long::GetOptions( "input_gff=s" => \$input_gff,
                          "output_gff=s" => \$output_gff,
                          "new_ids_list=s" => \$new_ids_list,
                          "help" => \$help) || die;

if (defined($help)) {
    print $usage;
    exit;
}

if (!defined($input_gff) ||  !(-e $input_gff) ||
    !defined($output_gff) ||  (-e $output_gff) ||
    !defined($new_ids_list) ||  !(-e $new_ids_list)) {
    die $usage;
}

my %old_to_new_ids;

open IN, $new_ids_list or die "\nUnable to open $new_ids_list for reading.\n\n";

while (my $line = <IN>) {

    chomp $line;
    my ($new_id, $old_id) = split "\t", $line;

    $old_to_new_ids{$old_id} = $new_id;

}

open OLD_GFF, $input_gff or die "\nUnable to open $input_gff for writing.\n\n";
open NEW_GFF, ">$output_gff" or die "\nUnable to open $output_gff for reading.\n\n";

while (my $line = <OLD_GFF>) {
    chomp $line;

    if ($line =~ /^#/) {
	print NEW_GFF "$line\n";
	next;
    }

    my @elems = split "\t", $line;
    if ($elems[1] eq 'maker' && $elems[2] eq 'gene' && $elems[8] =~ /ID=(.+);Name/) {
	my $id = $1;
	if (exists($old_to_new_ids{$id})) {
	    $elems[8] =~ s/$id/$old_to_new_ids{$id}/g;
	    my $new_line = join "\t", @elems;
	    print NEW_GFF "$new_line\n";
	}
	else {
	    print "$line\n";
	    die "\nUnable to find $id in old-to-new hash for gene.\n\n";
	}
    }
    elsif ($elems[1] eq 'maker' && $elems[2] eq 'mRNA' && $elems[8] =~ /ID=(.+);Parent/) {
	my $id = $1;

	my $base_id;
	my $isoform_num;
	if ($id =~ /(.+)-mRNA-(\d+)$/) {
	    $base_id = $1;
	    $isoform_num = $2;
	}
	else {
	    die "\nUnable to parse the gene id:  $id\n\n";
	}

	if (exists($old_to_new_ids{$base_id})) {
	    my $new_id = $old_to_new_ids{$base_id} . "." . $isoform_num;
	    $elems[8] =~ s/Parent=$base_id/Parent=$old_to_new_ids{$base_id}/;
	    $elems[8] =~ s/Name=$id/Name=$new_id/;
	    $elems[8] =~ s/ID=$id/ID=$new_id/;
	    my $new_line = join "\t", @elems;
	    print NEW_GFF "$new_line\n";
	}
	else {
	    print "$line\n";
	    die "\nUnable to find $id in old-to-new hash for mRNA.\n\n";
	}
    }
    elsif ($elems[1] eq 'maker' && $elems[2] eq 'CDS' && $elems[8] =~ /ID=(.+):cds;Parent/) {
	my $id = $1;

	my $base_id;
	my $isoform_num;
	if ($id =~ /(.+)-mRNA-(\d+)$/) {
	    $base_id = $1;
	    $isoform_num = $2;
	}
	else {
	    die "\nUnable to parse the CDS id:  $id\n\n";
	}

	if (exists($old_to_new_ids{$base_id})) {
	    my $new_id = $old_to_new_ids{$base_id} . "." . $isoform_num;
	    $elems[8] =~ s/Parent=$id/Parent=$new_id/;
	    $elems[8] =~ s/ID=$id/ID=$new_id/;
	    my $new_line = join "\t", @elems;
	    print NEW_GFF "$new_line\n";
	}
	else {
	    print "$line\n";
	    die "\nUnable to find $id in old-to-new hash for CDS.\n\n";
	}
    }
    elsif ($elems[1] eq 'maker' && $elems[2] eq 'exon' && $elems[8] =~ /ID=(.+):\d+;Parent/) {

	print "exon\n$line\n";
	
	my $id = $1;

	my $base_id;
	my $isoform_num;
	if ($id =~ /(.+)-mRNA-(\d+)/) {
	    $base_id = $1;
	    $isoform_num = $2;
	    print "$base_id\t$isoform_num\n";
	}
	else {
	    die "\nUnable to parse the exon id:  $id\n\n";
	}
#Cfim_AS236_tig06	maker	exon	2427134	2429137	.	-	.	ID=augustus_masked-Cfim_AS236_tig06-processed-gene-24.19-mRNA-1:exon:1664;Parent=augustus_masked-Cfim_AS236_tig06-processed-gene-24.19-mRNA-1
#augustus_masked-Cfim_AS236_tig06-processed-gene-24.19	1
#id augustus_masked-Cfim_AS236_tig06-processed-gene-24.19-mRNA-1:exon
#new id Cfim_AS236_t06g05070.1
#NEW:  Cfim_AS236_tig06	maker	exon	2427134	2429137	.	-	.	ID=Cfim_AS236_t06g05070.1:1664;Parent=augustus_masked-Cfim_AS236_tig06-processed-gene-24.19-mRNA-1


	if (exists($old_to_new_ids{$base_id})) {
	    my $new_id = $old_to_new_ids{$base_id} . "." . $isoform_num;

	    my $orig_exon_id = $id;
	    $orig_exon_id =~ s/\:exon//;

	    print "id $id\n";
	    print "orig_exon_id $orig_exon_id\n";
	    print "new id $new_id\n";
	    
	    $elems[8] =~ s/Parent=$orig_exon_id/Parent=$new_id/;
	    $elems[8] =~ s/ID=$id/ID=$new_id/;
	    my $new_line = join "\t", @elems;

	    print "NEW:  $new_line\n";
	    
	    print NEW_GFF "$new_line\n";
	}
	else {
	    print "$line\n";
	    die "\nUnable to find $id in old-to-new hash for exon.\n\n";
	}
    }
    elsif ($elems[1] eq 'maker' && $elems[2] =~ /prime_UTR/ && $elems[8] =~ /ID=(.+):\w+;Parent/) {
	my $id = $1;

	my $base_id;
	my $isoform_num;
	if ($id =~ /(.+)-mRNA-(\d+)$/) {
	    $base_id = $1;
	    $isoform_num = $2;
	}
	else {
	    die "\nUnable to parse the UTR id:  $id\n\n";
	}

	if (exists($old_to_new_ids{$base_id})) {
	    my $new_id = $old_to_new_ids{$base_id} . "." . $isoform_num;
	    $elems[8] =~ s/Parent=$id/Parent=$new_id/;
	    $elems[8] =~ s/ID=$id/ID=$new_id/;
	    my $new_line = join "\t", @elems;
	    print NEW_GFF "$new_line\n";
	}
	else {
	    print "$line\n";
	    die "\nUnable to find $id in old-to-new hash for UTR.\n\n";
	}
    }

    elsif ($elems[1] eq 'maker') {
	print "Should we consider this feature?\n";
	print "$line\n";
	print NEW_GFF "$line\n";
    }
    else {
	print NEW_GFF "$line\n";
    }
}

exit;

