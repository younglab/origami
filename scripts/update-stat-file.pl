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
    
    $stats{$opt} = $val;
  }
  
  close($sh);
}

switch($param) {
  case ["origamiver","date","cmdline","mode"] { $stats{$param} = join(" ",@args); }
  case ["firstreadfiles","secondreadfiles"] { $stats{$param} = join(",",@args); }
  case ["initialpets","procpets","bridgewithlinker","bridgewithoutlinker","macs1peaks","macs2peaks"] { 
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
      case "AB-linker" { }
    }

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
    print S "$key: $stats{$key}\n";
  }
  
  close(S);
}





