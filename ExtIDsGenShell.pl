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

my ($in, $dir, $prefix, $number, $name, $ref, $sn, $help);

GetOptions
(
	"i|in=s"=>\$in,
	"p|pf=s"=>\$prefix,
	"d|dir=s"=>\$dir,
	"m|nm=s"=>\$name,
	"n|num=i"=>\$number,
	"r|ref=s"=>\$ref,
	"sn=i"=>\$sn,
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
		-r <string> <the absolute path of genome reference> (dlft: /home/luna/Desktop/database/homo_minimap2/hsa.fa)
		-sn <int> <decide whether to run the later analysis according to the split of the bam file, ${red}recommend set to 1${end}> (default: 1)
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
	$ref = "/home/luna/Desktop/database/homo_minimap2/hsa.fa";
}

die "${red}$ref !exists$end\n" if(!(-e "$ref"));

my $refdir = dirname(File::Spec->rel2abs( $ref ));
$refdir=~s/^\s+|\s+$//g;

if(!(defined $sn))
{
	$sn = 1;
}

my $samtools;
if( -e "/home/luna/Desktop/Software/samtools/samtools/samtools" )
{
	$samtools = "/home/luna/Desktop/Software/samtools/samtools/samtools";
	$samtools=~s/^\s+|\s+$//g;
}
elsif(`which samtools`)
{
	$samtools = `which samtools`;
	chomp $samtools;
	$samtools=~s/^\s+|\s+$//g;
}
else
{
	die "${red}not samtools in your enviroment$end\n";
}

my $sambam;
if(-e "/home/luna/Desktop/Software/sambamba/build/sambamba")
{
	$sambam = "/home/luna/Desktop/Software/sambamba/build/sambamba";
	$sambam=~s/^\s+|\s+$//g;
}
elsif(`which sambamba`)
{
	$sambam = `which sambamba`;
	chomp $sambam;
	$sambam=~s/^\s+|\s+$//g;
}
else
{
	die "${red}not sambamba in your enviroment$end\n";
}

my $header = "${in}.header.sam";

if($in =~ /\.bam$/){
	open IN, "$samtools view -@ 3 $in |" || die $!;
	`$samtools view -H $in > $header`;
}
elsif($in =~ /\.sam\.gz$/){
	open IN, "gzip -dc $in |" || die $!;
	`gzip -dc $in | head -n 10000 | grep "^@" > $header`;
}
else{
	open IN, "$in |" || die $!;
	`head -n 10000 $in | grep "^@" > $header`;
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
my $sam = "NA";

#$samtools = "/home/luna/Desktop/Software/samtools/samtools-1.8/samtools";

while(my $line = <IN>){
	chomp $line;

	next unless($line !~ /^\@/);
	my @F = split /\t/, $line;

	if($idnum % $number == 0 && !exists $id{$F[0]}){
		close BAM;
		my $part = int($idnum/$number) + 1;
		$split = $part;
		open OB, "> $dir/$prefix.$name.part.$part" || die $!;
		# `$sambam index -t 20 $bam &` if($bam);

		`$samtools view -h -Sb --reference $ref --threads 10 -o $bam $sam && rm -rf $sam &` if(-e $sam && $bam);
		$bam = "$dir/$prefix.$name.part.$part.bam";
		$sam = "$dir/$prefix.$name.part.$part.sam";
		`cp $header $sam`;
		open BAM, " >> $sam" || die $!;
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

`$samtools view -h -Sb --reference $ref --threads 10 -o $bam $sam && rm -rf $sam &` if(-e $sam && $bam);

open SH, "> $dir/$name.runFind.sh" || die $!;
#print SH "#!/bin/bash\n";

print "${red}$name have $idnum read, the ids were split to $split part, complete this step$end\n";

for my $i(1..$split){
	$bam = "$dir/$prefix.$name.part.$i.bam";
#	`sambamba index -t 20 $bam`;
	print SH "#!/bin/bash\n";
	print SH "nohup perl $find -i $bam -o $dir/${name}.part.${i}.chimera -r $refdir &> $dir/$name.part.${i}.log &\nsleep 5m\n" if($sn);
	print SH "nohup perl $find -i $in -id $dir/$prefix.$name.part.$i -o $dir/${name}.part.${i}.chimera -r /home/luna/Desktop/database/homo_minimap2 &> $dir/$name.part.${i}.log &\nsleep 3m\n" if(! $sn);
}

close SH;

die "${red}$dir/$name.runFind.sh is empty, Please Check${end}\n" if(!(-s "$dir/$name.runFind.sh"));

if(-s "$dir/$name.runFind.sh"){
	print "${red}NOTICE: chmod +x $dir/$name.runFind.sh${end}\n";
	`chmod +x $dir/$name.runFind.sh`;
	print "${red}NOTICE: Want to find chimeras, just run the command \" $dir/$name.runFind.sh \" for starting${end}\n";
}
