#!/usr/bin/perl

use strict;

if( scalar(@ARGV) < 3 ) {
  die "binned-genome.pl <genome anno> <bin size> <min reads> <output prefix>";
}

my ($anno,$binsize,$minreads,$prefix) = @ARGV;

unless( -e $anno ) {
  die "binned-genome.pl: need genome annotation";
}

unless($binsize > 0 ) {
  die "binned-genome.pl: need positive bin size";
}

unless($minreads > 0 ) {
  die "binned-genome.pl: need positive minimum number of reads";
}

open(A,"<",$anno) or die "Cannot read $anno: $!";

my %chromosomes;

while(<A>) {
  chomp;
  my ($chr,$size) = split/\t/;
  
  $chromosomes{$chr} = $size;
}

close(A);

open(T,">","$prefix/tmp.bed") or die "Cannot open tmp.bed for writing: $!";

for my $key (keys(%chromosomes)) {
  my $bpsize = $chromosomes{$key} + $binsize;
  for(my $i = 0; $i <= $bpsize; $i += $binsize) {
    my $e = $i+$binsize-1;
    print T "$key\t$i\t$e\n";
  }
}

close(T);

### need to update this to handle different versions of bedtools
`bedtools coverage -counts -abam mapped_reads.bam -b $prefix/tmp.bed > $prefix/coverage.bed`;

open(F,"<","$prefix/coverage.bed") or die;
open(O,">","$prefix/coverage-peaks.bed") or die;

while(<F>) {
  chomp;
  
  my ($chr,$s,$e,$c) = split /\t/;
  
  print "$chr\t$s\t$e\n" unless $c < $minreads;
}

close(O);
close(F);