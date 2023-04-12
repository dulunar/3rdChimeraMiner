#!/usr/bin/perl -w
#########################################################################
# FileName: LapLengthGC.pl
# Version: 3c471015-3ab5-4178-93c3-372472b54b3b
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Tue 05 Jan 2021 04:27:48 PM CST
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;

my $red = "\033[0;31m";
my $end = "\033[0m";

my (@ins, $out, $samp, $min_lap, $max_lap, $help);
GetOptions
(
	"i|in=s"=>\@ins,
	"o|out=s"=>\$out,
	"mx=i"=>\$max_lap,
	"mi=i"=>\$min_lap,
	"help|?"=>\$help,
	"n|nm=s"=>\$samp,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file1> -i <file2> ...  -i <fileN> <in.input files>
		-o <string> <out.output file>
		-n <string> <the name of this sample>
		-mx <int> <the max overlap length> ${red}dflt = 62$end
		-mi <int> <the min overlap length> ${red}dflt = 0$end
INFO

die $usage if ($help || !@ins || !$out || !$samp);

if(! defined $min_lap)
{
	$min_lap = 0;
}

if(! defined $max_lap)
{
	$max_lap = 62;
}

open OUT,"> $out" || die $!;

my %number;
my %statgc;

my $tag;
my $tag1;
foreach my $file(@ins)
{
	chomp $file;
	open IN,"$file" || die $!;
	while(my $a = <IN>)
	{
		my $b=<IN>;
		my $c=<IN>;
		my $line=<IN>;
		my $d=<IN>;
		my $e=<IN>;
		chomp $line;
		my @F = split /\t/, $b;
		my @R = split /\t/, $c;
		
		$tag = ($F[1] ne $R[1]) ? "I" : "D";
		
		if(!$tag1)
		{
			$tag1 = $tag;
		}

		last if($tag1 ne $tag && $file =~ /inverted|direct/);

		$tag1 = $tag;

		my ($seq1, $nt, $dist) = (split /\t/,$line)[0,1,2];
		$seq1 =~ s/overlap:\s+//g;
		chomp $nt;
		$nt =~ s/\s+nt//g;

		my ($seq) = (split /;/, $seq1)[-1];
		chomp $seq;
		$seq =~ s/\s+nt//g;

		my $gc = $seq =~ tr/GCgc//;

		if(exists $number{$nt})
		{
			$number{$nt} ++ ;
			$statgc{$nt} += $gc;
		}
		else
		{
			$number{$nt} = 1;
			$statgc{$nt} = $gc;
		}
	}
	close IN;
	if($file !~ /inverted|direct/)
	{
		$tag = "T";
	}
}

print OUT "LapLength-$tag\tNumber-${samp}\tAll_Bases\tGC_Bases\t%GC-${samp}\n";
foreach my $nt ($min_lap..$max_lap)
{
	$number{$nt} = 0 if(!exists $number{$nt});
	$statgc{$nt} = 0 if(!exists $statgc{$nt});
	my $all_nt = $number{$nt}*$nt;
	my $per = 0;
	$per = ($all_nt > 0) ? $statgc{$nt}/$all_nt : 0;
	print OUT "${nt}\t$number{$nt}\t$all_nt\t$statgc{$nt}\t$per\n";
}
close OUT;
