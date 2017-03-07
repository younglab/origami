#!/usr/bin/perl

use strict;

die "Need to pass a bamfile" unless (@ARGV) >= 1;

my ($bamfile) = @ARGV;

open(S,"samtools view -H $bamfile| " ) or die "Cannot extract BAM header for $bamfile: $!";

while(<S>) {
  chomp;
  
  next unless /\@SQ.*SN:(\w+).*LN:(\w+)/;
  
  print "$1\t$2\n";
}

close(S);