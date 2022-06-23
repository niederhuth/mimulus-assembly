#! /opt/perl/bin/perl

# rename_gff_with_new_identifier.pl

# 15 October 2014

# Kevin Childs

# Given a maker gff3 file with awkward maker identifiers, make a new gff3 file with
# new identifiers.  A file with the new-to-old identifiers is also required.

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

    $old_id =~ s/-mRNA-1//;
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
    elsif ($elems[1] eq 'maker' && $elems[2] eq 'mRNA' && $elems[8] =~ /Parent=([^;]+)/) {
	my $id = $1;
	if (exists($old_to_new_ids{$id})) {
	    $elems[8] =~ s/Parent=$id/Parent=$old_to_new_ids{$id}/;
	    $elems[8] =~ s/Name=$id/Name=$old_to_new_ids{$id}/;
	    $elems[8] =~ s/ID=$id/ID=$old_to_new_ids{$id}/;
	    my $new_line = join "\t", @elems;
	    print NEW_GFF "$new_line\n";
	}
	else {
	    print "$line\n";
	    die "\nUnable to find $id in old-to-new hash for mRNA.\n\n";
	}
    }
    elsif ($elems[1] eq 'maker') {
	print NEW_GFF "$line\n";
    }
    else {
	print NEW_GFF "$line\n";
    }
}

exit;

