#!/usr/bin/perl

use strict; 

my ($results) = @ARGV;
my %anchors;

open(F,"<",$results) or die "Cannot open file $results: $!";

<F>; #ignore header

while(<F>) {
  chomp;
  
  my ($chr1,$s1,$e1,$chr2,$s2,$e2) = split /,/;
  
  my $val1 = "$chr1\t$s1\t$e1";
  my $val2 = "$chr2\t$s2\t$e2";
  $anchors{$val1} = 1;
  $anchors{$val2} = 1;
}

for my $key (keys(%anchors)) {
  print "$key\n";
}

close(F);