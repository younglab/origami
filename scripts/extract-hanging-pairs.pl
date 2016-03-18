#!/usr/bin/perl

use strict;

if( $#ARGV < 2 ) {
	print "file both one neither\n";
	exit 1;
}

my $outboth = $ARGV[1];
my $outone = $ARGV[2];
my $outneither = $ARGV[3];

my $header= `samtools view -H $ARGV[0]`;

open(F,"samtools view $ARGV[0] |") or die "Cannot open: $!";
open(OB,">$outboth.sam") or die "cannot open: $!";
open(OO1,">$outone" . "_1.fq") or die "cannot open: $!";
open(OO2,">$outone" . "_2.fq") or die "cannot open: $!";
open(ON1,">$outneither" ."_1.fq") or die "cannot open: $!";
open(ON2,">$outneither" ."_2.fq") or die "cannot open: $!";


print OB $header;

until(eof(F)) {
	my $line1 = <F>;
	my $line2 = <F>;

	chomp $line1;
	chomp $line2;

	my @a1 = split /\t/, $line1;
	my @a2 = split /\t/, $line2;

	if( $a1[2] ne "*" && $a2[2] ne "*" ) {
		print OB "$line1\n";
		print OB "$line2\n";
	} elsif( $a1[2] eq "*" && $a2[2] ne "*") {
		print OO1 "\@$a1[0]\n$a1[9]\n+$a1[0]\n$a1[10]\n";
                print OO2 "\@$a2[0]\n$a2[9]\n+$a2[0]\n$a2[10]\n";
	} elsif( $a2[2] eq "*" && $a1[2] ne "*") {
                print OO1 "\@$a2[0]\n$a2[9]\n+$a2[0]\n$a2[10]\n";
                print OO2 "\@$a1[0]\n$a1[9]\n+$a1[0]\n$a1[10]\n";
	} else {
                print ON1 "\@$a1[0]\n$a1[9]\n+$a1[0]\n$a1[10]\n";
                print ON2 "\@$a2[0]\n$a2[9]\n+$a2[0]\n$a2[10]\n";
	}
}

`samtools view -Sb $outboth.sam > $outboth.bam`;

close(OB);
close(OO1);
close(OO2);
close(ON1);
close(ON2);
close(F);
