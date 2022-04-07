#!/usr/bin/perl -w
#########################################################################
# FileName: FindChimeras.pl
# Version: 919f20ee-6aeb-4722-9e52-d55526be0b45
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Wed 20 Oct 2021 08:45:00 PM CST
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;
my $red = "\033[0;31m";	my $end = "\033[0m";
my ($in,$out,$refdir,$id,$percent_mis,$help);
GetOptions
(
	"i|in=s"=>\$in,
	"o|out=s"=>\$out,
	"p|ps=f"=>\$percent_mis,
	"r|rd=s"=>\$refdir,
	"d|id=s"=>\$id,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file> <in.input file>
		-o <string> <out.output file>
		-id <string> <the list of QNAMEs splitted for this part> (Optical)
		${red}-p <float> <the allowed mismatch of overlap sequence> (Default: 0.2)$end
		${red}-r <directory for storing the reference and each chromosome> <default: /home/luna/Desktop/database/homo_bwa>$end
INFO

die $usage if ($help || !$in || !$out);

if(!(defined $percent_mis)){
	$percent_mis = 0.2;
}

if(!(defined $refdir)){
	$refdir = "/home/luna/Desktop/database/homo_bwa";
}

my $samtools;
if( -e "/home/luna/Desktop/Software/samtools/samtools/samtools" ){
	$samtools = "/home/luna/Desktop/Software/samtools/samtools/samtools";
}
elsif(`which samtools`){
	$samtools = `which samtools`;
	chomp $samtools;
	$samtools=~s/^\s+|\s+$//g;
}
else{
	die "not samtools in /home/luna/Desktop/Software/samtools/samtools/samtools\n";
}

my %mis = ();
for my $len(0..61){
	my $ok = $len * $percent_mis;
	my $o;
	if($ok=~/^\d+$/){
		if($ok == 0){
			$o = 1;
		}
		else{
			$o = $ok;
		}
	}
	else{
		$o = int($ok) + 1;
	}
	$mis{$len} = $o;
}

die "$id do not exists\n" if(defined $id && !(-e "$id"));

if($in =~ /\.bam$/){
	if(! defined($id)){
		open IN, "$samtools view $in | " || die $!;
	}
	else{
		open IN, "$samtools view -@ 5 -N $id $in |" || die $!;
	}
}
elsif($in =~ /\.sam$/){
	if(! defined($id)){
		open IN, "$in" || die $!;
	}
	else{
		open IN, "fgrep -w -f $id $in |" || die $!;
	}
}
elsif($in =~ /\.sam\.gz$/){
	if(! defined($id)){
		open IN, "gzip -dc $in | " || die $!;
	}
	else{
		open IN, "gzip -dc $in | fgrep -w -f $id | " || die $!;
	}
}

open OUT,"> $out" || die $!;

my %read;
my $detail;
my @sa = ();

while(my $line = <IN>){
	chomp $line;
	my @F = split /\t/, $line;
	
	#next unless($F[2] !~ /M|MT/);

	my $flag = reverse( sprintf("%b",$F[1]) );
	my @tmp = split //, $flag;
	next unless(!$tmp[8] || $tmp[8] == 0);

	if(!exists $read{$F[0]}){
		$read{$F[0]} = 1;
		@sa = ();
	}
	else{
		$read{$F[0]} ++;
	}

	next unless($read{$F[0]} <= 2);

	if($read{$F[0]} == 1){
		$detail = $line;
	}
	else{
		my @FF = split /\t/, $detail;

		my $str = &strand_code($FF[1]);
		my $read_seq = ($str eq "+") ? ($FF[9]) : (&seq_trans($FF[9]));
		my $read_seq_rc = ($str eq "+") ? (&seq_trans($FF[9])) : ($FF[9]);

		my $SA1;
		for my $i(@FF){
			if($i =~ /SA:Z:/){
				($SA1) = $i =~ /SA:Z:(.*)/;
			}
		}

		my $SA2;
		for my $i(@F){
			if($i =~ /SA:Z:/){
				($SA2) = $i =~ /SA:Z:(.*)/;
			}
		}

		my @sa = split /;/, $SA1;
		my @sa2 = split /;/, $SA2;
		unshift(@sa, $sa2[0]);

		my %sat = ();
		my @bmaps = ();
		my $i = -1;
		foreach my $saa (@sa){
			$i ++;
			my ($chr, $pos_plus, $str_seg, $cigar, $mapq, $num_mis) = split /,/, $saa;
			my $map;
			my $bmap;
			my $indel;
			($cigar, $bmap, $map, $indel) = &callengthcigar($str_seg, $cigar);
			push @bmaps, $bmap;
			$sat{$i} = "$chr,$pos_plus,$str_seg,$cigar,$bmap,$map,$indel";
		}

		my ($b, $j) = &bmapsort(\@bmaps);
		my @sbmap = @$b;
		my @tags = @$j;
		
		my $read_length = length $read_seq;

		print OUT "$FF[0]\t$read_length\t$read_seq\n";
		my $sn = 0;

		my $segment_ref_sequence;
		my $reverse_extend_forward;
		my $reverse_extend_end;
		my $tag_lap = "D";
		my $cut_forward_seq;
		my $cut_end_seq;
		my $ref_reads_sequence;
		my ($overlap_forward, $overlap_forward_length) = ("", 0);
		my ($overlap_end, $overlap_end_length) = ("", 0);
		my ($overlap_seq, $overlap_seq_length) = ("", 0);
	
		$sn ++;
		my ($end_plus, $sub_start, $seg_start, $seg_end, $seg_seq) = ("", "", "", "", "", "");
		my ($chr, $pos_plus, $str_seg, $cigar, $bmap, $map, $indel) = split /,/, $sat{$tags[0]};
		$end_plus = $pos_plus + $map - 1;
		if($str_seg eq "+"){
			$sub_start = ($bmap == 0) ? 0 : $bmap;
			$seg_start = ($bmap == 0) ? 1 : $bmap + 1;
			my $seg_len = $map + $indel;
			$seg_end = $seg_start + $seg_len - 1;
			$seg_seq = substr($read_seq, $sub_start, $seg_len);
			print OUT "Seg_$sn\tOverlap\tDistance\t$chr\t$str_seg\t$pos_plus\t$end_plus\t$seg_start\t$seg_end\t$seg_seq\n";
		}
		else{
			$sub_start = ($bmap == 0) ? 0 : $bmap;
			$seg_start = ($bmap == 0) ? 1 : $bmap + 1;
			my $seg_len = $map + $indel;
			$seg_end = $seg_start + $seg_len - 1;
			$seg_seq = substr($read_seq, $sub_start, $seg_len);
			print OUT "Seg_$sn\tOverlap\tDistance\t$chr\t$str_seg\t$end_plus\t$pos_plus\t$seg_start\t$seg_end\t$seg_seq\n";
		}

		for my $i(1..$#tags){
			$sn ++;
			my ($sub_start1, $seg_start1, $seg_end1, $seg_seq1) = (0, 0, 0, "");

			my ($chr1, $pos_plus1, $str_seg1, $cigar1, $bmap1, $map1, $indel1) = split /,/, $sat{$tags[$i]};
			my $end_plus1 = $pos_plus1 + $map1 - 1;
			$sub_start1 = ($bmap1 == 0) ? 0 : $bmap1;
			$seg_start1 = ($bmap1 == 0) ? 1 : $bmap1 + 1;
			my $seg_len1 = $map1 + $indel1;
			$seg_end1 = $seg_start1 + $seg_len1 - 1;
			$seg_seq1 = substr($read_seq, $sub_start1, $seg_len1);

			my $tmp_reads_referecne_end = $pos_plus + $map - 1;
			$ref_reads_sequence = &ref_reads_seq($chr, $pos_plus, $tmp_reads_referecne_end);

			if($seg_start1 > $seg_end){
				if($str_seg eq "+"){	
					my $tmp_start = $end_plus + 1;
					my $tmp_end = $end_plus + 31;
					$segment_ref_sequence = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$reverse_extend_forward = &length_cut_2(uc($segment_ref_sequence));
					#print "$chr, $str_seg, $tmp_start, $tmp_end";
					$tmp_start = $end_plus - 30;
					$tmp_end = $end_plus;
					$cut_forward_seq = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$cut_forward_seq = &length_cut_1(uc($cut_forward_seq));

					#my $cut_seg = substr($seg_seq, -31);
					#print "\t$chr, $str_seg, $tmp_start, $tmp_end";

					if($str_seg1 eq "+"){
						$tmp_start = $pos_plus1;
						$tmp_end = $pos_plus1 + 30;
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(uc($cut_end_seq));
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";

						$tmp_start = $pos_plus1 - 31;
						$tmp_end = $pos_plus1 - 1;
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$reverse_extend_end = &length_cut_1(uc($segment_ref_sequence));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length = length $overlap_seq;
						my $distance = $pos_plus1 - $end_plus;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$pos_plus1\t$end_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
					else{
						$tmp_start = $end_plus1 - 30;
						$tmp_end = $end_plus1;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(&seq_trans(uc($cut_end_seq)));

						$tmp_end = $end_plus1 + 31;
						$tmp_start = $end_plus1 + 1;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$reverse_extend_end = &length_cut_1(&seq_trans(uc($segment_ref_sequence)));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length =
						length $overlap_seq;
						my $distance = $pos_plus1 - $end_plus - $overlap_seq_length;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$end_plus1\t$pos_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
				}
				else{
					my $tmp_end = $pos_plus - 1;
					my $tmp_start = $pos_plus - 31;
					#print "$chr, $str_seg, $tmp_start, $tmp_end";
					$segment_ref_sequence = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$reverse_extend_forward = &length_cut_2(&seq_trans(uc($segment_ref_sequence)));

					$tmp_end = $pos_plus + 30;
					$tmp_start = $pos_plus;
					#print "\t$chr, $str_seg, $tmp_start, $tmp_end";
					$cut_forward_seq = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$cut_forward_seq = &length_cut_2(&seq_trans(uc($cut_forward_seq)));

					#my $cut_seg = substr($seg_seq, -31);
					if($str_seg1 eq "+"){
						$tmp_start = $pos_plus1;
						$tmp_end = $pos_plus1 + 30;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(uc($cut_end_seq));

						$tmp_start = $pos_plus1 - 31;
						$tmp_end = $pos_plus1 - 1;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$reverse_extend_end = &length_cut_1(uc($segment_ref_sequence));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length = length $overlap_seq;
						my $distance = $pos_plus1 - $end_plus - $overlap_seq_length;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$pos_plus1\t$end_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
					else{
						$tmp_end = $end_plus1;
						$tmp_start = $end_plus1 - 30;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(&seq_trans(uc($cut_end_seq)));

						$tmp_end = $end_plus1 + 31;
						$tmp_start = $end_plus1 + 1;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$reverse_extend_end = &length_cut_1(&seq_trans(uc($segment_ref_sequence)));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length = length $overlap_seq;
						my $distance = $pos_plus - $end_plus1;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$end_plus1\t$pos_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
				}
			}
			## when seg_start1 < seg_end, it needs some change for overlap
			else{
				my $overlap_seg_len = $seg_end - $seg_start1 + 1;
				my $overlap_seg_start = $seg_start1 - 1;

				if($str_seg eq "+"){	
					my $tmp_start = $end_plus + 1;
					my $tmp_end = $end_plus + 31;
					$segment_ref_sequence = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$reverse_extend_forward = &length_cut_2(uc($segment_ref_sequence));
					#print "$chr, $str_seg, $tmp_start, $tmp_end";
					$tmp_start = $end_plus - 30;
					$tmp_end = $end_plus;
					$cut_forward_seq = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$cut_forward_seq = &length_cut_1(uc($cut_forward_seq));

					#my $cut_seg = substr($seg_seq, -31);
					#print "\t$chr, $str_seg, $tmp_start, $tmp_end";

					if($str_seg1 eq "+"){
						$tmp_start = $pos_plus1 + $overlap_seg_len;
						$tmp_end = $pos_plus1 + $overlap_seg_len + 30;
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(uc($cut_end_seq));
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";

						$tmp_start = $pos_plus1 + $overlap_seg_len - 31;
						$tmp_end = $pos_plus1 + $overlap_seg_len - 1;
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$reverse_extend_end = &length_cut_1(uc($segment_ref_sequence));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length = length $overlap_seq;
						my $distance = $pos_plus1 - $end_plus;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$pos_plus1\t$end_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
					else{
						$tmp_start = $end_plus1 - $overlap_seg_len - 30;
						$tmp_end = $end_plus1 - $overlap_seg_len;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(&seq_trans(uc($cut_end_seq)));

						$tmp_end = $end_plus1 - $overlap_seg_len + 31;
						$tmp_start = $end_plus1 - $overlap_seg_len + 1;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$reverse_extend_end = &length_cut_1(&seq_trans(uc($segment_ref_sequence)));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length =
						length $overlap_seq;
						my $distance = $pos_plus1 - $end_plus - $overlap_seq_length;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$end_plus1\t$pos_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
				}
				else{
					my $tmp_end = $pos_plus - 1;
					my $tmp_start = $pos_plus - 31;
					#print "$chr, $str_seg, $tmp_start, $tmp_end";
					$segment_ref_sequence = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$reverse_extend_forward = &length_cut_2(&seq_trans(uc($segment_ref_sequence)));

					$tmp_end = $pos_plus + 30;
					$tmp_start = $pos_plus;
					#print "\t$chr, $str_seg, $tmp_start, $tmp_end";
					$cut_forward_seq = &ref_reads_seq($chr, $tmp_start, $tmp_end);
					$cut_forward_seq = &length_cut_2(&seq_trans(uc($cut_forward_seq)));

					#my $cut_seg = substr($seg_seq, -31);

					if($str_seg1 eq "+"){
						$tmp_start = $pos_plus1 + $overlap_seg_len;
						$tmp_end = $pos_plus1 + $overlap_seg_len + 30;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(uc($cut_end_seq));

						$tmp_start = $pos_plus1 + $overlap_seg_len - 31;
						$tmp_end = $pos_plus1 + $overlap_seg_len - 1;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$reverse_extend_end = &length_cut_1(uc($segment_ref_sequence));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length = length $overlap_seq;
						my $distance = $pos_plus1 - $end_plus - $overlap_seq_length;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$pos_plus1\t$end_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
					else{
						$tmp_end = $end_plus1 - $overlap_seg_len;
						$tmp_start = $end_plus1 - $overlap_seg_len - 30;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end";
						$cut_end_seq = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$cut_end_seq = &length_cut_2(&seq_trans(uc($cut_end_seq)));

						$tmp_end = $end_plus1 - $overlap_seg_len + 31;
						$tmp_start = $end_plus1 - $overlap_seg_len + 1;
						#print "\t$chr1, $str_seg1, $tmp_start, $tmp_end\n";
						$segment_ref_sequence = &ref_reads_seq($chr1, $tmp_start, $tmp_end);
						$reverse_extend_end = &length_cut_1(&seq_trans(uc($segment_ref_sequence)));

						#my $cut_seg1 = substr($seg_seq1, 0, 31);
						#print "seg1	$cut_seg1\t";
						($overlap_end, $overlap_end_length) = &compare_sequence_end($cut_end_seq, $reverse_extend_forward);
						#print "seg	$cut_seg\t";
						($overlap_forward, $overlap_forward_length) = &compare_sequence_forward($cut_forward_seq, $reverse_extend_end);

						$overlap_seq = $overlap_forward.$overlap_end;

						$overlap_seq_length = length $overlap_seq;
						my $distance = $pos_plus - $end_plus1;
						my $over_fromat = "$overlap_seq:${overlap_seq_length}nt";
						print OUT "Seg_$sn\t$over_fromat\t$distance\t$chr1\t$str_seg1\t$end_plus1\t$pos_plus1\t$seg_start1\t$seg_end1\t$seg_seq1\n";
					}
				}

			}

			($end_plus, $sub_start, $seg_start, $seg_end, $seg_seq) = ($end_plus1, $sub_start1, $seg_start1, $seg_end1, $seg_seq1);
			($chr, $pos_plus, $str_seg, $cigar, $bmap, $map, $indel) = ($chr1, $pos_plus1, $str_seg1, $cigar1, $bmap1, $map1, $indel1);

		}
	}
}
close IN;
close OUT;

sub compare_sequence_forward {
	my ($cut_forward_seq, $reverse_extend_end) = @_;
	#print "$cut_forward_seq\t$reverse_extend_end\n";
	my $overlap = "";
	my $length_overlap = 0;
	my @cut_read_seq = split //,$cut_forward_seq;
	my @extend_ref_seq = split //,$reverse_extend_end;
	my $fault = 0;
	my $len = 0;
	my $num = $#cut_read_seq;
	while($num >= 0) {
		$overlap = $cut_read_seq[$num].$overlap;
		if($cut_read_seq[$num] ne $extend_ref_seq[$num]){
			$fault ++;
		}
		$len = ($overlap ne "") ? length($overlap) : 0;
		if($fault > $mis{$len}){
			if($cut_read_seq[$num+1] ne $extend_ref_seq[$num+1]){
				$overlap = substr($overlap, 2);
			}
			elsif($cut_read_seq[$num+1] eq $extend_ref_seq[$num+1]){
				$overlap = substr($overlap, 1);	
			}
		}
		last if($fault > $mis{$len});
		$num --;
	}

	$length_overlap = length $overlap;
	return($overlap, $length_overlap);
}

sub compare_sequence_end{
	my ($cut_end_seq, $reverse_extend_forward) = @_;
	#print "$cut_end_seq\t$reverse_extend_forward\n";
	my $overlap = "";
	my $length_overlap = 0;
	my @cut_read_seq = split //,$cut_end_seq;
	my @extend_ref_seq = split //,$reverse_extend_forward;
	my $fault = 0;
	my $num = 0;
	my $len = 0;
	while($num <= $#cut_read_seq){
		$overlap .= $cut_read_seq[$num];
		if($cut_read_seq[$num] ne $extend_ref_seq[$num]){
			$fault ++;
		}
		$len = ($overlap) ? length($overlap) : 0;
		if($fault > $mis{$len}){
			if($cut_read_seq[$num-1] eq $extend_ref_seq[$num-1]){
				$overlap = substr($overlap, 0, -1);	
			}
			elsif($cut_read_seq[$num-1] ne $extend_ref_seq[$num-1]){
				$overlap = substr($overlap, 0, -2);
			}
		}
		last if $fault > $mis{$len};
		$num ++;
	}
	$length_overlap = length $overlap;
	return($overlap, $length_overlap);
}

sub bmapsort {
	my ($bmaps) = @_;
	my @j;

	my @sort = sort { $a <=> $b } @$bmaps;
	for my $sbmap (@sort){
		my $i = -1;
		for my $bmap (@$bmaps){
			$i ++;
			push @j, $i if( $sbmap == $bmap && !( grep {$_ == $i} @j ) );
		}
	}
	return(\@sort, \@j);
}

sub callengthcigar {
	my ($str, $cigars) = @_;
	my @len = split /\D/, $cigars;
	my @type = split /\d+/, $cigars;
	my $cigar;
	if($str eq "-"){
		for(my $i=$#len; $i >= 0; $i --){
			$cigar .= "$len[$i]$type[$i+1]";
		}
	}
	else{
		$cigar = $cigars;
	}

	@len = split /\D/, $cigar;
	@type = split /\d+/, $cigar;
	my $map = 0;
	my $bmap = 0;
	my $indel = 0;
	for(my $i=0; $i <= $#len; $i ++){
		$map += $len[$i] if($type[$i+1] =~ /M/);
		$indel += $len[$i] if($type[$i+1] =~ /I/);
	}

	for(my $i=0; $i <= $#len; $i ++){
		last if($type[$i+1] =~ /M|I/);
		$bmap += $len[$i] if($type[$i+1] =~ /S/);
	}

	return ($cigar, $bmap, $map, $indel);
}

sub strand_code {
	my $flag = $_[0]; #数组是@_，但是单个参数是$_[0]开始；
	chomp $flag;
	$flag = reverse( sprintf("%b",$flag) );
	my @tmp = split //, $flag;
	my $strand = (!$tmp[4] || $tmp[4] == 0) ? ("+") : ("-");
	return $strand;
}

sub seq_trans {
	my $seq = $_[0];
	chomp $seq;
	my $reverse_seq;
	$seq = uc($seq);
	$seq=~tr/ATCG/TAGC/;
	$reverse_seq = reverse $seq;
	return $reverse_seq;
}

sub seq_comp {
	my $seq = $_[0];
	chomp $seq;
	my $comp_seq;
	$seq = uc($seq);
	$seq=~tr/ATCG/TAGC/;
	$comp_seq = $seq;
	return $comp_seq;
}

sub ref_reads_seq{
	my ($chr,$start,$end) = @_;
	die "$refdir/$chr.fa did exists\n" if(! -e "$refdir/$chr.fa");
	`$samtools faidx $refdir/$chr.fa` if(! -e "$refdir/$chr.fa.fai");
	my $fai = `cat $refdir/$chr.fa.fai`;
	chomp $fai;
	my $length_chr = (split /\s+/, $fai)[1];
	$end = ($length_chr > $end) ? $end : $length_chr;
	$start = ($length_chr > $start) ? $start : $length_chr;
	my $ref_reads_sequence = `$samtools faidx $refdir/$chr.fa $chr:$start-$end`;
	my @ALL = split /\n/,$ref_reads_sequence; 
	$ref_reads_sequence=""; 
	for(my $i = 1; $i<=$#ALL; $i++){
		$ref_reads_sequence .= $ALL[$i];
	}
	$ref_reads_sequence = uc($ref_reads_sequence);
	return $ref_reads_sequence;
}

sub length_cut_1{
	my $seq_cut = $_[0];
	if(length($seq_cut) < 31){
		my $tmp_seq = "N"x(31 - length($seq_cut));
		$seq_cut = $tmp_seq.$seq_cut;
	}
	$seq_cut = uc($seq_cut);
	return $seq_cut;
}

sub length_cut_2{
	my $seq_cut = $_[0];
	if(length($seq_cut) < 31){
		my $tmp_seq = "N"x(31 - length($seq_cut));
		$seq_cut = $seq_cut.$tmp_seq;
	}
	$seq_cut = uc($seq_cut);
	return $seq_cut;
}

sub fuzzy_matching {
	my ($original_pattern, $mismatches_allowed) = @_;
	$mismatches_allowed >= 0 || die "Number of mismatches must be greater than or equal to zero\n";
	my $new_pattern = build_approximate($original_pattern, $mismatches_allowed);
	return($new_pattern);
}

sub build_approximate {
	my ($pattern, $mismatches_allowed) = @_;
	if ($mismatches_allowed == 0) {
		return $pattern;
	}
	elsif (length($pattern) <= $mismatches_allowed) {
		$pattern =~ tr/ACTG/./;
		return $pattern;
	}
	else {
		my ($first, $rest) = $pattern =~ /^(.)(.*)/;
		my $after_match = build_approximate($rest, $mismatches_allowed);
		if ($first =~ /[ACGT]/) {
			my $after_miss = build_approximate($rest, $mismatches_allowed - 1);
			return "(?:$first$after_match|.$after_miss)";
		}
		else {
			return "$first$after_match";
		}
	}
}

