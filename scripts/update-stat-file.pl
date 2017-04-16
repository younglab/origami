#!/usr/bin/perl

use strict;
use Switch;
    use Fcntl qw(:flock);

die "Not enough arguments" unless scalar(@ARGV)>=3;

my ($statfile,$param,@args) = @ARGV;

my %stats;

my $stattypes = [
  ["Origami version","origamiver"],
  ["Date and time of run","date"],
  ["Command line arguments","cmdline"],
  ["First read files","firstreadfiles"],
  ["Second read files","secondreadfiles"],
  ["Initial number of PETs","initialpets"],
  ["Number of PETs after initial processing","procpets"],
  ["Linker-trimming mode","mode"],
  #["firstreadbridgetrimmed","Number of PETs with bridge linker"]
];

for my $ref (@{$stattypes}) {
  $stats{$ref->[1]} = [$ref->[0],""];
}

unless( ! -e $statfile ) {
  open(my $sh,"<","$statfile") or die;
  flock($sh,LOCK_EX) or die "Could not acquire lock on $statfile: $!";

  while(<$sh>) {
    chomp;
    my ($opt,$val) = split /: /;
    
    for my $ref (@{$stattypes}) {
      if($ref->[0] eq $opt) {
        $stats{$ref->[1]}->[1] = $val;
        last; ## right now, silently ignore anything it doesn't understand
      }
    }
  }
  
  flock($sh,LOCK_UN);
  close($sh);
}

switch($param) {
  case ["origamiver","date","cmdline","mode"] { $stats{$param}->[1] = join(" ",@args); }
  case ["firstreadfiles","secondreadfiles"] { $stats{$param}->[1] = join(",",@args); }
  case ["initialpets","procpets"] { 
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
      $stats{$param}->[1] = $nlines/4;
    }
}

open(O,">","$statfile") or die "Cannot write to $statfile: $!";

for my $ref (@{$stattypes}) {
  my $opt = $ref->[1];
  print O "$stats{$opt}->[0]: $stats{$opt}->[1]\n";
}

close(O);
