#!/usr/bin/perl -w
#########################################################################
# FileName: SplitChimericRead.pl
# Version: 0617b6bd-1b01-43a1-b7ae-4fff23d5c7fb
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Mon 25 Apr 2022 03:33:50 PM CST
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $idfile, $ic, $out, $help);

GetOptions
(
	"i|in=s"=>\$in,
	"f|if=s"=>\$idfile,
	"c|ic=s"=>\$ic,
	"out=s"=>\$out,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file> <in.input raw bam file>
		-f <file> <input id file>
		-c <file> <input chimera file>
		-o <string> <out.output fastq file>
INFO

die $usage if ($help || !$in || !$idfile || !$ic || !$out);

if($ic =~ /\.gz$/)
{
	open IC, "gzip -dc $ic |" || die $!;
}
else
{
	open IC, "$ic" || die $!;
}

open OID, "> $idfile" || die $!;

my %num = ();
while(my $L1 = <IC>)
{
	my $L2 = <IC>;
	my $L3 = <IC>;
	my $L4 = <IC>;
	my $L5 = <IC>;
	my $L6 = <IC>;

	chomp $L1;
	my ($id)=(split /\t/, $L1)[0]=~/(.*)\:\d+nt/;
	if(!exists $num{$id})
	{
		$num{$id} = 1;
		print OID "$id\n";
	}
	else
	{
		$num{$id} ++;
	}
}
close IC;
close OID;

open IN, "samtools view -N $idfile $in | " || die $!;
open OUT,"> $out" || die $!;

if($ic =~ /\.gz$/)
{
    open IC, "gzip -dc $ic |" || die $!;
}
else
{
    open IC, "$ic" || die $!;
}

my $cuta = 0;
my $cutb = 0;
my $last = 0;

my $tag = "!";
my %hash = ();
my $temp_id;

while(my $line = <IN>){
	chomp $line;
	my @F = split /\t/, $line;
	my $id = $F[0];
	my $seq = $F[9];
	my $qua = $F[10];

	my $len = length $seq;

	#my $num = ($ic =~ /\.gz$/) ? `zgrep -w $id $ic | wc -l` : `grep -w $id $ic | wc -l`;
	my $num = $num{$id};
	
	if($temp_id && $temp_id ne $id)
	{
		%hash = ();
	}

	my $start = 1;
	my $end = $len;
	for my $i(1..$num)
	{
		my $L1 = <IC>;
		my $L2 = <IC>;
		my $L3 = <IC>;
		my $L4 = <IC>;
		my $L5 = <IC>;
		my $L6 = <IC>;

		chomp $L5;
		my ($st1, $ed1, $st2, $ed2) = (split /\s+/, $L5)[1] =~ /(\d+)->(\d+):(\d+)->(\d+)/;
		my $Len1 = $ed1 - $st1 + 1;
		my $Len2 = $ed2 - $st2 + 1;

		if(!exists $hash{$id}{$st1}{$ed1})
		{
			$hash{$id}{$st1}{$ed1} = 1;
		}
		else
		{
			$hash{$id}{$st1}{$ed1} ++;
		}

		if(!exists $hash{$id}{$st2}{$ed2})
		{
			$hash{$id}{$st2}{$ed2} = 1;
		}
		else
		{
			$hash{$id}{$st2}{$ed2} ++;
		}

		chomp $L2;
		my ($seqF) = (split /\s+/, $L2)[4];
		my $lenF = length $seqF;
		print "${red}$id Seg$i Length is error${end}\n" if($lenF != $Len1);

		chomp $L3;
		my ($seqL) = (split /\s+/, $L3)[4];
		my $lenL = length $seqL;
		my $j = $i + 1;
		print "${red}$id Seg$j Length is error${end}\n" if($lenL != $Len2);

		my $quaF = $tag x $lenF;
		my $quaL = $tag x $lenL;

		if($start >= $st1)
		{
			print OUT "\@${id}.Seg${i}.$st1-$ed1\n$seqF\n+\n$quaF\n" if($hash{$id}{$st1}{$ed1} == 1);
			if($st2 > $ed1)
			{
				my $sub_len = $st2 - $ed1 - 1;
				if($sub_len >= 30)
				{
					my $qua_sub = $tag x $sub_len;
					my $sub = substr($seq, $ed1, $sub_len);
					my $stm = $ed1 + 1;
					my $edm = $st2 - 1;
					print OUT "\@${id}.CutA${cuta}.$stm-$edm\n$sub\n+\n$qua_sub\n";
					${cuta} ++;
				}
				print OUT "\@${id}.Seg${j}.$st2-$ed2\n$seqL\n+\n$quaL\n" if($hash{$id}{$st2}{$ed2} == 1);
			}

			else
			{
				print OUT "\@${id}.Seg${j}.$st2-$ed2\n$seqL\n+\n$quaL\n" if($hash{$id}{$st2}{$ed2} == 1);
			}
		}
		else
		{
			my $sub_len = $st1 - $start;
			if($sub_len >= 30)
			{
				my $qua_sub = $tag x $sub_len;
				my $sub = substr($seq, $start - 1, $sub_len);
				my $stm = $start;
				my $edm = $st1 - 1;
				print OUT "\@${id}.CutB${cutb}.$stm-$edm\n$sub\n+\n$qua_sub\n";
				${cutb} ++;
			}
			print OUT "\@${id}.Seg${i}.$st1-$ed1\n$seqF\n+\n$quaF\n" if($hash{$id}{$st1}{$ed1} == 1);
			if($st2 > $ed1)
			{
				my $sub_len = $st2 - $ed1 - 1;
				if($sub_len >= 30)
				{
					my $qua_sub = $tag x $sub_len;
					my $sub = substr($seq, $ed1, $sub_len);
					my $stm = $ed1 + 1;
					my $edm = $st2 - 1;
					print OUT "\@${id}.CutA${cuta}.$stm-$edm\n$sub\n+\n$qua_sub\n";
					${cuta} ++;
				}
				print OUT "\@${id}.Seg${j}.$st2-$ed2\n$seqL\n+\n$quaL\n" if($hash{$id}{$st2}{$ed2} == 1);
			}
			else
			{
				print OUT "\@${id}.Seg${j}.$st2-$ed2\n$seqL\n+\n$quaL\n" if($hash{$id}{$st2}{$ed2} == 1);
			}
		}
		$start = $ed2 + 1;
	}
	$temp_id = $id;

	if($end > $start)
	{
		my $sub_len = $end - $start + 1;
		next unless $sub_len >= 30;
		my $qua_sub = $tag x $sub_len;
		my $sub = substr($seq, $start - 1, $sub_len);
		my $stm = $start;
		my $edm = $end;
		print OUT "\@${id}.Last${last}.$stm-$edm\n$sub\n+\n$qua_sub\n";
		$last ++;
	}
}

close IN;
close OUT;
close IC;
