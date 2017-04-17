#!/usr/bin/perl

use strict;
use Switch;

die "Not enough arguments" unless scalar(@ARGV)>=3;

my ($statfile,$param,@args) = @ARGV;

my %stats;

unless( ! -e $statfile ) {
  open(my $sh,"<","$statfile") or die;

  while(<$sh>) {
    chomp;
    my ($opt,$val) = split /: /;
    
    if( $val =~ /^ARRAY/ ) {
      my @vals = split /\t/, $val;
      shift @vals;
      
      $stats{$opt} = \@vals;
    } else {
      $stats{$opt} = $val;
    }
  }
  
  close($sh);
}

switch($param) {
  case ["origamiver","date","cmdline","mode"] { $stats{$param} = join(" ",@args); }
  case ["firstreadfiles","secondreadfiles"] { $stats{$param} = join(",",@args); }
  case ["initialpets","procpets","bridgewithlinker","bridgewithoutlinker","macs1peaks","macs2peaks",
  "ablinkeraa","ablinkerbb","ablinkerab","ablinkernottrimmed"] { 
      my $nlines = 0;
      if($args[0] eq "-" ) {
        while(<STDIN>) {
          print;
          $nlines++;
        }
      } else {
        for my $f (@args) {
          open(F,"<","$f") or die;
          $nlines++ while(<F>);
          close(F);
        }
      }
      $stats{$param} = $nlines/4;
    }
  case ["bowtiefirstread","bowtiesecondread"] {

    my ($proc,$reported,$failed,$rep) = (0,0,0,0);
    
    for my $f (@args) {
      open(B,"<","$f") or die "Cannot read bowtie summary file: $!";
      
      my @data = <B>;
      chomp @data;
      my ($p,$a,$f,$r) = @data;
      
      $p =~ /(\d+)/;
      $proc += $1;
      $a =~ /(\d+)/;
      $reported += $1;
      $f =~ /(\d+)/;
      $failed += $1;
      $r =~ /(\d+)/;
      $rep += $1;
      
      close(B);
    }
    
    $stats{$param} = [$proc,$reported,$failed,$rep];
  }
  case "final" {
    open(O,">","$args[0]") or die "Cannot write final report to $args[0]: $!";
    
    print O <<END;
Run report for origami alignment
================================

Origami version: $stats{"origamiver"}
Date and time of start of run: $stats{"date"}
Command arguments used: $stats{"cmdline"}
First read files: $stats{"firstreadfiles"}
Second read files: $stats{"secondreadfiles"}

Initial number of PETs: $stats{"initialpets"}
Number of PETs after pre-processing: $stats{"procpets"}

Linker-trimming run mode: $stats{"mode"}
END
  
    switch($stats{"mode"}) {
      case "bridge-linker" { 
        print O<<END;
Number of PETs with bridge-linker: $stats{"bridgewithlinker"}
Number of PETs without bridge-linker: $stats{"bridgewithoutlinker"}
END
        }
      case "AB-linker" { 
        print O<<END;
Number of PETs with A linker on both ends: $stats{"ablinkeraa"}
Number of PETs with B linker on both ends: $stats{"ablinkerbb"}
END
        print O "Number of PETs with same linker sequence on both ends (non-chimeric): " . ($stats{"ablinkeraa"}+$stats{"ablinkerbb"}) . "\n";
        print O "Number of PETs with different linker sequences on both ends (chimeric): $stats{'ablinkerab'}\n";
        print O "Number of PETs that lacked any linker sequence: $stats{'ablinkernottrimmed'}\n";
        }
    }
    
    print O<<END;

Number of first end reads in bowtie aligment: $stats{"bowtiefirstread"}->[0]
Number of first end reads mapping to the reference genome: $stats{"bowtiefirstread"}->[1]
Number of first end reads failing to map to the reference genome: $stats{"bowtiefirstread"}->[2]
Number of first end reads discarded because they repetitively mapped to the reference genome: $stats{"bowtiefirstread"}->[3]


Number of second end reads in bowtie aligment: $stats{"bowtiesecondread"}->[0]
Number of second end reads mapping to the reference genome: $stats{"bowtiesecondread"}->[1]
Number of second end reads failing to map to the reference genome: $stats{"bowtiesecondread"}->[2]
Number of second end reads discarded because they repetitively mapped to the reference genome: $stats{"bowtiesecondread"}->[3]
END
    print O<<END;
    
Number of MACS1 peaks called: $stats{"macs1peaks"}
Number of MACS2 peaks called: $stats{"macs2peaks"}
END
    close(O);
  }
}

unless($param eq "final") {
  open(S,">","$statfile") or die "Cannot update stat file $statfile: $!";
  
  for my $key (keys(%stats)) {
    unless(ref($stats{$key}) eq "ARRAY") {
      print S "$key: $stats{$key}\n";
    } else {
      print S "$key: ARRAY\t" . join("\t",@{$stats{$key}}) . "\n";
    }
  }
  
  close(S);
}





