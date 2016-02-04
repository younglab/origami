#!/usr/bin/perl

### 10, 11, 12

use strict;

my %c = {};

open(F,$ARGV[0]) or die "Cannot open $ARGV[0]: $!";

while(<F>) {
	chomp;
	my @l1 = split /\t/;
	$_ = <F>;
	chomp;
	my @l2 = split /\t/;

	my ($c1,$p11,$p12) = @l1[10,11,12];
	my ($c2,$p21,$p22) = @l2[10,11,12];

	next if( $c1 =~ /chrM/ or $c2 =~ /chrM/);

	my $out = "";

#	print STDERR "$c1,$p11,$p12 $c2,$p21,$p22\n";

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

for my $k (sort keys(%c)) {
	next if $k =~ /HASH/; ## weird error...
	my @arr = split /:/, $k;
	print "$_\t" for (@arr);
	print "$c{$k}\n";
}

close(F);
