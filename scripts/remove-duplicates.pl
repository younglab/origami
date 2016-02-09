#!/usr/bin/perl

use strict;

open(F,"samtools view -h $ARGV[0]|") or die "Cannot open $ARGV[0]: $!";

my %h;
my $n = 0;
my $ex = 0;
my $line1;
my $line2;

while(<F>) {
  $line1 = $_;
  last unless /^@/;
  print;
}

chomp $line1;
$line2 = <F>;
chomp $line2;
while(!eof(F)) {
  my @a1 = split /\t/, $line1;
  my @a2 = split /\t/, $line2;
  $n++;
  # concatenate the reference name and position as a key
  my $s = "$a1[2]:$a1[3]:$a2[2]:$a2[3]";
  if(!exists($h{$s})) {
    $h{$s} = 1;
    print "$line1\n$line2\n";
  } else {
    $ex++;
    #print STDERR "Excluded $a1[0]:$a2[0]\n";
  }
  $line1 = <F>;
  $line2 = <F>;
  chomp $line1;
  chomp $line2;
}

print STDERR "Processed $n PETs, excluded $ex PETs because of duplication\n";

close(F);