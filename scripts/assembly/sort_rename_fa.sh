#!/bin/bash

#Ugly script to sort and rename fasta files

cat $1 |\
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1rn |\
cut -f 2- |\
tr "\t" "\n" |\
awk -v i=00001 '/^>/{print ">contig" ++i, $0; next}{print}' 

