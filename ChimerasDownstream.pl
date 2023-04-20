#!/usr/bin/perl -w
#########################################################################
# FileName: ChimerasDownstream.pl
# Version: c90c1b8b-a9b1-47a7-be52-93938028ec91
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Thu 09 Dec 2021 03:37:05 PM CST
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;

my $red = "\033[0;31m";
my $end = "\033[0m";

my $usage = "
Usage: perl $0

${red}For Extracting, transformating, and getinfo of chimeras.$end

	-d <string> <the directory of all chimeras file>
	-od <string> <the directory for saving the results> <dflt = \$dir/chimeras_analysis_lap(the min overlap length)$end>
	-n <string> <the name of this sample>
	-L <INT> <min length of segment> <dflt = 50>
	-p <INT> <the min overlap length> <dflt = 3>
	-s <INT> <the smallest distance of two segment> <dflt = 25>
	-b <INT> <the biggest distance of two segment> <dflt = 10000>

Example: ${red}perl $0 -d /home/luna/work/ThirdChimera/0907data/5e.2 -n 5e.2 -L 50 -s 25 -b 10000$end\n";

my ($dir, $name, $minLen, $help, $molpl, $smd, $bgd, $outdir);

GetOptions(
	'd=s'		=>	\$dir,
	'od=s'		=>	\$outdir,
	'help|?'	=>	\$help,
	'L=i'		=>	\$minLen,
	'n=s'		=>	\$name,
	'p=i'		=>	\$molpl,
	's=i'		=>	\$smd,
	'b=i'		=>	\$bgd,
);

die "$usage\n" if ($help || !$dir || !$name);

if(!(defined $minLen)){
	$minLen = 50;
}

if(!(defined $molpl)){
	$molpl = 2;
}

if(!(defined $smd)){
	$smd = 25;
}

if(!(defined $bgd)){
	$bgd = 10000;
}

if(!(defined $outdir)){
	if($molpl == 3){
		$outdir = "$dir/chimeras_analysis_lap${molpl}";
	}
	else{
		$outdir = "$dir/chimeras_analysis_lap${molpl}";
	}
}

my $split = 0;
if(-s "$dir/$name.part.1.chimera"){
	$split = `ls $dir/$name.part.*.chimera | wc -l`;
}

die "${red}$dir/$name.part.1.chimera !exists$end\n" if($split == 0);

`mkdir -p $outdir` if(!(-d "$outdir"));

my $out = "$outdir/$name.error.txt";
my $normal = "$outdir/$name.normal.txt";

open OUT, "> $out" || die $!;
open NM, "> $normal" || die $!;

open OT,"> $outdir/$name.TransFormat.txt" || die $!;

open OD,"> $outdir/$name.direct" || die $!;
open OV,"> $outdir/$name.inverted" || die $!;

open OI,"> $outdir/$name.chimera.id" ||die $!;

my ($invert, $direct) = (0) x 2;
my ($La, $Lb, $Lc, $Ld) = (0) x 4;

my $tag;
my $in;
my $len_total = 0;
my $chimera_site = 0;
my $read = 0;

my %hash = ();

for my $i (1..$split){
	if(-s "$dir/$name.part.$i.chimera"){
		$in = "$dir/$name.part.$i.chimera";
	}
	elsif(-s "$dir/$name.part.$i.chimera.gz"){
		$in = "$dir/$name.part.$i.chimera.gz";
	}
	else{
		die "no file in this directory $dir\n";
	}

	next unless (defined $in);

	if($in =~ /\.gz$/){
		open IN,"gzip -dc $in |" || die $!;	
	}
	else{
		open IN, "$in" || die $!;
	}

	my $line = <IN>;
	chomp $line;

	my ($id, $length, $seq) = split /\t/, $line;
	my ($id1, $length1, $seq1) = split /\t/, $line;

	$len_total += $length;

	$line = <IN>;
	chomp $line;
	my ($seg, $lap, $dis, $chr, $str, $gstart, $gend, $rstart, $rend, $segseq) = split /\t/, $line;
	my ($seg1, $lap1, $dis1, $chr1, $str1, $gstart1, $gend1, $rstart1, $rend1, $segseq1) = split /\t/, $line;

	my @seg = ();

	my $tran = 0;
	my $tran_dc = 0;
	my $tran_iv = 0;
	my $tran_dr = 0;

	while (my $line = <IN>) {
		chomp $line;

		if($line !~ /^Seg/){
			($seg1, $lap1, $dis1, $chr1, $str1, $gstart1, $gend1, $rstart1, $rend1, $segseq1) = ($seg, $lap, $dis, $chr, $str, $gstart, $gend, $rstart, $rend, $segseq);
			($id, $length, $seq) = split /\t/, $line;
			$len_total += $length;
		}
		else{
			($seg, $lap, $dis, $chr, $str, $gstart, $gend, $rstart, $rend, $segseq) = split /\t/, $line;
			if($id eq $id1){
				$tran ++;
				my ($obas, $fseqp, $ref_f, $lseqp, $ref_l, $seqp, $lenp) = $lap =~ /(.*):(.*);(.*):(.*);(.*):(.*):(.*)nt/;

				#$dis = $gend1 - $gstart - $lenp;
				my $trans = "$id:${length}nt\t$seg1->$seg\n$chr1\t$str1\t$gstart1\t -> $segseq1 -> \t$gend1\n$chr\t$str\t$gstart\t -> $segseq  -> \t$gend\noverlap: $obas;$fseqp:$ref_f;$lseqp:$ref_l;$seqp\t$lenp nt\tdistance: $dis nt\nstart ${rstart1}->${rend1}:${rstart}->${rend} end\n\n";

				my $len1 = length($segseq1);
				my $len = length($segseq);
				
				print OT "$trans";

				if($chr eq $chr1){
					if ($len < $minLen || $len1 < $minLen){
						print OUT "$trans";
					}
					elsif( abs($dis) > $bgd || abs($dis) < $smd ){
						print OUT "$trans";
					}
					elsif( abs($lenp) < $molpl || $lenp eq "" ){
						print OUT "$trans";
					}
					#	elsif($rstart - $rend1 > 50)
					#{
					#	print OUT "$trans";
					#}
					else{
						print NM "$trans";

						if(!exists $hash{$id}){
							$hash{$id} = 1;
							$read ++;
							print OI "$id\n";
						}
						else{
							$hash{$id} ++;
						}

						$chimera_site ++;
						if ($str1 eq "+" && $str eq "-"){
							$La ++;
							$tag = "invert";
							$invert ++;
							print OV "$trans";
						}
						elsif ($str1 eq "-" && $str eq "+"){
							$Lb ++;
							$tag = "invert";
							$invert ++;
							print OV "$trans";
						}
						elsif ($str1 eq "+" && $str eq "+"){
							$Lc ++;
							$tag = "direct";
							$direct ++;
							print OD "$trans";
						}
						elsif ($str1 eq "-" && $str eq "-"){
							$Ld ++;
							$tag = "direct";
							$direct ++;
							print OD "$trans";
						}

					}
				}
				else{
					print OUT "$trans";
				}
				($seg1, $lap1, $dis1, $chr1, $str1, $gstart1, $gend1, $rstart1, $rend1, $segseq1) = ($seg, $lap, $dis, $chr, $str, $gstart, $gend, $rstart, $rend, $segseq);
			}
			else{
				($id1, $length1, $seq1) = ($id, $length, $seq);
				($seg1, $lap1, $dis1, $chr1, $str1, $gstart1, $gend1, $rstart1, $rend1, $segseq1) = ($seg, $lap, $dis, $chr, $str, $gstart, $gend, $rstart, $rend, $segseq);
			}
		}
	}
	close IN;
}

close OUT;
close OT;
close OV;
close OD;
close OI;
close NM;

open OO,"> $outdir/$name.chimera.id.times" ||die $!;

foreach my $id(keys %hash){
	print OO "$id\t$hash{$id}\n";
}
close OO;

print "Title\t+-\t-+\t++\t--\tInvertChimera\tDirectChimera\tReadsNumber\tChimericSite\tTotalLengthofReads\n";
print "$name\t$La\t$Lb\t$Lc\t$Ld\t$invert\t$direct\t$read\t$chimera_site\t$len_total\n";
