#!/usr/bin/perl

use strict;

my $file = $ARGV[0];
my $dc = $ARGV[1];

open(F,"$file") or die "Cannot open $file! $!";

<F>; #skip header line

while(<F>) {
  chomp;
  my @arr = split /,/;
  
  next unless $arr[0] eq $arr[3]; # Remove interchromosomal interactions
  next if $arr[0] =~ /chrM/; # the browser doesn't like chrM
  
  print "$arr[0]:$arr[1]-$arr[2]\t$arr[3]:$arr[4]-$arr[5]\t$arr[5+$dc]\n";
}

close(F);