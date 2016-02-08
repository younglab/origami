#!/usr/bin/perl

### 10, 11, 12

use strict;

my %c;

open(F,$ARGV[0]) or die "Cannot open $ARGV[0]: $!";

my $line1 = <F>;
my $line2 = <F>;

while(!eof(F)) {
	chomp $line1;
	my @l1 = split /\t/, $line1;
	chomp $line2;
	my @l2 = split /\t/, $line2;

  if(!($l1[0] eq $l2[0] && $l1[1] == $l2[1] && $l1[2] == $l2[2] &&
      $l1[3] eq $l2[3] && $l1[4] == $l2[4] && $l1[5] == $l2[5])) {
      ## somehow we got disoriented, toss a line a move forward
      
      $line1 = $line2;
      $line2 = <F>;
      next;
  }

	my ($c1,$p11,$p12) = @l1[10,11,12];
	my ($c2,$p21,$p22) = @l2[10,11,12];
	
	
	### exclude chrM (although it would be nice to generalize this exclusion to arbitrary chromosomes)
	unless( $c1 =~ /chrM/ or $c2 =~ /chrM/) {

	  my $out = "";


	  if( $c1 gt $c2 ) {
		  $out = "$c2:$p21:$p22:$c1:$p11:$p12";
	  }
	  elsif( $c1 eq $c2 ) {
		  if( $p11 > $p21 ) {
	       $out = "$c2:$p21:$p22:$c1:$p11:$p12";
		  } else {
	                $out = "$c1:$p11:$p12:$c2:$p21:$p22";
		  }
	  } else {
                $out = "$c1:$p11:$p12:$c2:$p21:$p22";

	  }
	 
	  $c{$out}++;
	}
	$line1 = <F>;
	$line2 = <F>;
}

for my $k (sort keys(%c)) {
	next if $k =~ /HASH/; ## weird error...
	my @arr = split /:/, $k;
	print "$_\t" for (@arr);
	print "$c{$k}\n";
}

close(F);
