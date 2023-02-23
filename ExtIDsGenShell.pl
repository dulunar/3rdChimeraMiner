#!/usr/bin/perl -w
#########################################################################
# FileName: ExtIDsGenShell.pl
# Version: a9250741-d72f-4239-89fb-dbb5324566f4
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Mon 06 Dec 2021 10:09:09 PM CST
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $red = "\033[0;31m";
my $end = "\033[0m";

my $homepath = dirname(File::Spec->rel2abs( $0 ));
$homepath=~s/^\s+|\s+$//g;
my $find = "$homepath/FindChimeras.latest.pl";
die "${red}$find do not exists${end}\n" if(!(-e "$find"));

my ($in, $dir, $prefix, $number, $name, $ref, $help);

GetOptions
(
	"i|in=s"=>\$in,
	"p|pf=s"=>\$prefix,
	"d|dir=s"=>\$dir,
	"m|nm=s"=>\$name,
	"n|num=i"=>\$number,
	"r|ref=s"=>\$ref,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		${red}-i <file> <in.input file mem.bam|mem.sam, must not sorted by chromosome>$end
		-p <string> <the prefix of out.output file> (Default: id)
		-d <string> <the root directory for storing the file of \$name fold, it contains id list and chimera>
		-m <string> <the name of this sample>
		${red}-n <INT> <the maximal reads id number in one list> (Defaults: 100000)$end
		-r <string> <the absolute path of genome reference> (dlft: /home/luna/Desktop/database/homo_bwa/hsa.fa)
INFO

die $usage if ($help || !$in || !$dir || !$name);

my $indir = dirname(File::Spec->rel2abs( $in ));
$indir=~s/^\s+|\s+$//g;

my $file = basename $in;
$file=~s/^\s+|\s+$//g;

$in = "$indir/$file";

die "${red}$in !exists$end\n" if(!(-e "$in"));

if(!(defined $number))
{
	$number = 100000;
}

if(!(defined $prefix))
{
	$prefix = "id";
}

if(!(defined $ref))
{
	$ref = "/home/luna/Desktop/database/homo_bwa/hsa.fa";
}

if($in =~ /\.bam$/){
	open IN, "samtools view -@ 3 $in |" || die $!;
}
elsif($in =~ /\.sam\.gz$/){
	open IN, "gzip -dc $in |" || die $!;
}
else{
	open IN, "$in |" || die $!;
}

if(!(-e "$dir/$name")){
	`mkdir -p $dir/$name`;
}
$dir = "$dir/$name";

open OUT,"> $dir/$prefix.$name.txt" || die $!;

my $split = 0;

my %id = ();

my $idnum = 0;

my $bam;

while(my $line = <IN>){
	chomp $line;

	next unless($line !~ /^\@/);
	my @F = split /\t/, $line;

	if($idnum % $number == 0 && !exists $id{$F[0]}){
		my $part = int($idnum/$number) + 1;
		$split = $part;
		open OB, "> $dir/$prefix.$name.part.$part" || die $!;
		`sambamba index -t 20 $bam &` if($bam);
		$bam = "$dir/$prefix.$name.part.$part.bam";
		open BAM, " | samtools view -bS --reference $ref --threads 4 -o $bam - " || die $!;
	}

	if(!exists $id{$F[0]}){
		$id{$F[0]} = 1;
		$idnum ++;
		print OUT "$F[0]\n";
		print OB "$F[0]\n";
		print BAM "$line\n";
	}
	else{
		$id{$F[0]} ++;
		print BAM "$line\n";
	}
}
close OUT;
close IN;

open SH, "> $dir/$name.runFind.sh" || die $!;
print SH "#!/bin/bash\n";
print "${red}$name have $idnum read, the ids were split to $split part$end\n";
for my $i(1..$split){
	$bam = "$dir/$prefix.$name.part.$i.bam";
#	`sambamba index -t 20 $bam`;
	print SH "nohup perl $find -i $dir/$bam -o $dir/${name}.part.${i}.chimera -p 0.2 -r /home/luna/Desktop/database/homo_bwa &> $dir/$name.part.${i}.log &\nsleep 3m\n" if($ref);
	print SH "nohup perl $find -i $in -id $dir/$prefix.$name.part.$i -o $dir/${name}.part.${i}.chimera -p 0.2 -r /home/luna/Desktop/database/homo_bwa &> $dir/$name.part.${i}.log &\nsleep 3m\n" if(! $ref);
}

close SH;

die "${red}$dir/$name.runFind.sh is empty, Please Check${end}\n" if(!(-s "$dir/$name.runFind.sh"));

if(-s "$dir/$name.runFind.sh"){
	print "${red}NOTICE: chmod +x $dir/$name.runFind.sh${end}\n";
	`chmod +x $dir/$name.runFind.sh`;
	print "${red}NOTICE: Want start. just run the command \" $dir/$name.runFind.sh \"${end}\n";
}
