#!/usr/bin/perl

use strict;

my ($resultsfile,$anchorfile) = @ARGV;
my %anchors;
my %anchorscore;

open(R,"<",$resultsfile) or die "Cannot open $resultsfile: $!";
open(A,"<",$anchorfile) or die "Cannot open $anchorfile: $!";

while(<A>) {
  chomp;
  
  my ($achr,$as,$ae,$pchr,$ps,$pe,undef,$pscore,undef) = split /\t/;
  
  my $anchor = "$achr\t$as\t$ae";
  my $peak = "$pchr,$ps,$pe";
  
  $pscore = 0 if !defined($pscore);
  
  ## for now, take first hit as the "peak"
  if(!exists($anchors{$anchor}) || $anchorscore{$anchor} < $pscore) {
    $anchors{$anchor} = $peak;
    $anchorscore{$anchor} = $pscore;
  }
}

my $header = <R>;
print $header;

while(<R>) {
  chomp;
  
  my ($chr1,$s1,$e1,$chr2,$s2,$e2,@val) = split /,/;

  my $r1 = $anchors{"$chr1\t$s1\t$e1"};
  my $r2 = $anchors{"$chr2\t$s2\t$e2"};

  print "$r1,$r2";
  print ",$_" for(@val);
  print "\n";
}

close(R);
close(A);