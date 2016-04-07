#!/usr/bin/perl

use strict;
use Switch;

my ($mode,$base,$ablinker,$outfile) = @ARGV;

sub dirmode {
  my $b = shift;
  my $o = shift;
  
  my $leftb = 0;
  my $leftt = 0;
  my $leftc = 0;
  my $leftu = 0;
  my $leftd = 0;
  my $c = 0;

  
  opendir(my $dirh, $b);
  
  while(readdir($dirh)) {
    next unless /leftreads/;
    
    open(F,"<","$b/$_") or die "cannot open $_: $!";
    
    $c = 0;
    $c++ while(<F>);
    
    $leftb += $c/4;
    
    close(F);
  }
  closedir($dirh);

  
  open(F,"<","$b/left_kept.fq") or die "Cannot open $b/left_kept.fq: $1";
  $leftt++ while(<F>);
  close(F);
  
  open(F,"<","$b/left_untrimmed.fq") or die "Cannot open $b/left_untrimmed.fq: $!";
  $leftu++ while(<F>);
  close(F);
  
  $leftt /= 4;
  $leftu /= 4;
  
  if( $ablinker =~ /yes/ ) {
    open(F,"<","$b/left_chimeric.fq") or die "Cannot open $b/left_chimeric.fq: $!";
    $leftc++ while(<F>);
    close(F);
    
    $leftc /= 4;
  } 
  
  $leftd = $leftb - $leftt - $leftc - $leftu;

  open(O,">",$outfile) or die "cannot open $outfile";
  
  print O "Trimming statistics:\n";
  print O "Total read pairs: $leftb\n";
  if( $ablinker =~ /yes/ ) {
    print O "Trimmed AA/BB read pairs: $leftt\n";
    print O "Trimmed AB/BA read pairs: $leftc\n";
  } else {
    print O "Trimmed read pairs: $leftt\n";
  }
  print O "Untrimmed read pairs: $leftu\n";
  print O "Pairs discarded because of trimmed read less than set minimum length: $leftd\n";
  
  close(O);
}

sub filemode {
  
}

switch($mode) {
  case "dir" { dirmode($base,$outfile) }
  case "file" { filemode($base,$outfile) }
}