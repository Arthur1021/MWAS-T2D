#!/usr/bin/perl -w

use strict;
use warnings;

my $cpgList = $ARGV[0];
my $dir = $ARGV[1] || die "Usage: CpG.list.tab output.folder\n";
my %cpg; # $cpg{name}{chr/pos}

open (my $fh, $cpgList) || die "Can't open $cpgList: $!\n";
while (<$fh>) {
    chomp;
    my @d = split(/\t/, $_);
    $cpg{$d[0]}{chr} = $d[1];
    # $cpg{$d[0]}{pos} = $d[2];
}
close $fh;

opendir (D, $dir) || die "Can't open dir $dir: $!\n";
my @f = grep { /\.chr/ && /\.tab/ } readdir D;
closedir D;

foreach my $cpg (sort keys %cpg) {
    my $chr = $cpg{$cpg}{chr};
    # my $dir = 'extracted.for.metal';

    my $fc = 0;
    my $str = "SCHEME\tSTDERR\n";
    foreach my $f (@f) {
	if ($f =~ /$cpg/) {
	$fc++;
	# print "find $f for $cpg\n";
	# $str .= "\nMARKER\tsnps\nALLELE\tsnps.coded snps.noncoded\nEFFECT\tbeta\nSTDERR\tstd.err\nPVAL\tpvalue\n";
	$str .= "\nMARKER\tsnps\nALLELE\tsnps.coded snps.noncoded\nEFFECT\tbeta\nSTDERR\tstd.err\n"; # pvalue is not needed when SCHEME STDERR is used 
	$str .= "PROCESS\t$f\n";
	}
    }
    # $str .= "\nOUTFILE\tmeta.$cpg.out\n"; # output file name does not work; have to rename metal outputs for each cpg or make subfolders
    $str .= "\nOUTFILE\tmetal.$cpg. .tbl\nANALYZE\n";
    
    if ($fc == 1) {
	my $outFile = $dir . '/single.metal.' . $cpg . '.txt';
	open (my $ofh, ">$outFile") || die "Can't open $outFile: $!\n";
	print $ofh $str;
	close $ofh;
    } else {
	# reset
	print STDERR "$cpg has $fc .tab data files\n";
    }
}
