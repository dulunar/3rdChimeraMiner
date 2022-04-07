#!/usr/bin/perl -w
#########################################################################
# FileName: QC_normalized.pl
# Version: 28578195-be8e-44ac-b6a8-aae5a6138ef5
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Thu 24 Mar 2022 12:20:48 PM CST
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;
my ($in, $out, $val, $help);
GetOptions
(
	"i|in=s"=>\$in,
	"o|out=s"=>\$out,
	"q|qua=i"=>\$val,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file> <in.input fastq file>
		-o <string> <out.output fastq file>
		-q <INT> <the normalization quality value> (Default: 41)
INFO

die $usage if ($help || !$in || !$out);

if(!(defined $val)){
	$val = 41;
}

my $time = `date`;
print "Start Analysis: $time";

if($in =~ /\.gz$/) {
	open IN, "gzip -dc $in |" || die $!;
}
else{
	open IN, "$in" || die $!;
}

if($out =~ /\.gz$/) {
	open OUT,"| gzip > $out " || die $!;
}
else{
	open OUT,"| gzip > $out.gz " || die $!;
}


my %hash = ();
my $read = 0;
while(<IN>) {
	<IN>;
	<IN>;
	my $qua = <IN>;
	chomp $qua;
	$read ++;
	my @F = split //, $qua;
	my $i = 0;
	for my $r(@F) {
		$i ++;
		my $q = ord($r) - 33;
		if(!exists $hash{$i}) {
			$hash{$i} = $q;
		}
		else{
			$hash{$i} += $q;
		}
	}

	print "$read finished\n" if($read % 1000000 == 0);
}
close IN; 

$time = `date`;
print "Finish read fastq and Start cal mean value: $time";

print "ReadSite\tMeanQuality\n";
my %mean = ();
foreach my $i (sort {$a<=>$b} keys %hash) {
	my $mean = $hash{$i} / $read;
	$mean{$i} = $mean;
	print "$i\t$mean\n";
}

$time = `date`;
print "Finish cal mean value and Start modif fastq: $time";

$read = 0;
if($in =~ /\.gz$/) {
    open IN, "gzip -dc $in |" || die $!;
}
else{
    open IN, "$in" || die $!;
}

use POSIX;
while(my $a = <IN>) {
	my $b = <IN>; 
	my $c = <IN>;
	my $qua = <IN>;
	chomp $qua;
	my @F = split //, $qua;
	my $i = 0;
	my $line;
	$read ++;
	for my $r(@F) {
		$i ++; 
		my $q = ord($r) - 33;
		my $real = round($q * $val / $mean{$i}) + 33;
		my $rq = chr($real);
		$line .= $rq;
	}
	
	print OUT "$a$b$c$line\n";

	if($read % 1000000 == 0){
		$time = `date`;
		print "$read modified: $time";
	}
}
close IN;
close OUT;

$time = `date`;
print "Analysis Complete: $time";
