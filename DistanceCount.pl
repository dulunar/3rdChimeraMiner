#!/usr/bin/perl -w
#########################################################################
# FileName: DistanceCount.pl
# Version: 15756763-f752-4f09-bc64-d20e74d21ded
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Sun 10 Jan 2021 10:42:41 PM CST
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;

my ($in,$outdir,$out,$step_length,$help);
GetOptions
(
	"i|in=s"	=>	\$in,
	"d|od=s"	=>	\$outdir,
	"o|out=s"	=>	\$out,
	"l|sl=i"	=>	\$step_length,
	"help|?"	=>	\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file> <in.input file>
		-d <string> <output directory>
		-o <string> <the prefix of out.output file>
		-sl <int> <step length> dflt = 5
INFO

die $usage if ($help || !$in || !$out);
open IA, " $in " || die $!;

if(!(defined $outdir)){
	$outdir = dirname $in;
	chomp $outdir;
	$outdir=~s/^\s+|\s+$//g;
}
my $file = basename $in;
my $samp = (split /\./, $file)[0];

if(!(defined $step_length)){
	$step_length = 5;
}

my %hash=();
while(my $a=<IA>){
	my $b=<IA>;
	my $c=<IA>;
	my $d=<IA>;
	my $e=<IA>;
	<IA>;
	chomp $b;
	chomp $c;
	chomp $d;
	my @F = split /\t/, $d;
	my ($lap) = $F[0] =~ /overlap: (.*)/;
	my ($len) = $F[1] =~ /(\d+) nt/;
	my ($dil) = $F[2] =~ /distance: (.*) nt/;
	my $dis = int(abs($dil)/$step_length);
	if(exists $hash{$len}{$dis}){
		$hash{$len}{$dis} ++;
	}
	else{
		$hash{$len}{$dis} = 1;
	}
}
close IA;

my %out = ();
my $head = "distance";
my $max = int(10000/$step_length);
open OA,"> $outdir/$out.dis.$step_length.txt" || die $!;
foreach my $len(sort {$a<=>$b} keys %hash){
	$head .= "\t$len";
	foreach my $dis(0..$max){
		my $aa = (!exists $hash{$len}{$dis}) ? 0 : $hash{$len}{$dis};
		if(!exists $out{$dis}){
			$out{$dis} = $aa;
		}
		else{
			$out{$dis} .= "\t$aa";
		}
	}
}

print OA "$head\tAllLap\n";
foreach my $dis(sort {$a<=>$b} keys %out){
	chomp $out{$dis};
	my @F = split /\t/, $out{$dis};
	my $sum = 0;
	foreach my $i(@F){
		$sum += $i;
	}
	my $tmp = int($sum);
	print OA "$dis\t$out{$dis}\t$tmp\n";
}
close OA;
