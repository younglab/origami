#!/usr/bin/perl

use strict;
use Switch;

my ($mode,$base,$outfile) = @ARGV;

sub dirmode {
  my $b = shift;
  my $o = shift;
  
  my $leftb = 0;
  my $leftt = 0;
  my $leftu = 0;
  my $c = 0;

  
  opendir(my $dirh, $b);
  
  while(readdir($dirh)) {
    next unless /leftreads/;
    
    open(F,"<","$b/$_") or die "cannot open $_";
    
    $c = 0;
    $c++ while(<F>);
    
    $leftb += $c/4;
    
    close(F);
  }
  closedir($dirh);

  
  open(F,"<","$b/left_kept.fq") or die "Cannot open $b/left_kept.fq";
  $leftt++ while(<F>);
  
  close(F);
  
  open(F,"<","$b/left_untrimmed.fq") or die "Cannot open $b/left_kept.fq";
  $leftu++ while(<F>);
  close(F);

  open(O,">",$outfile) or die "cannot open $outfile";
  
  print O "Trimming statistics:\n";
  print O "Total read pairs: $leftb\n";
  print O "Trimmed read pairs: $leftt\n";
  print O "Untrimmed read pairs: $leftu\n";
  
  close(O);
}

sub filemode {
  
}

switch($mode) {
  case "dir" { dirmode($base,$outfile) }
  case "file" { filemode($base,$outfile) }
}