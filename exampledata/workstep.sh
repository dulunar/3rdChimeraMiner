#!/bin/bash
#########################################################################
# FileName: workstep.sh
# Version: 96f1bb65-c1b4-4857-84cc-016f62589e75
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Wed 12 Apr 2023 09:48:50 PM CST
#########################################################################

dir=`pwd` 
sdir=`dirname $dir`

ref=/home/luna/Desktop/database/homo_minimap2/hsa.fa
mmi=/home/luna/Desktop/database/homo_minimap2/hsa.HIFI.mmi
samtools=/home/luna/Desktop/Software/samtools/samtools/samtools
minimap2=/home/luna/Desktop/Software/minimap2/minimap2

# step 1, convert the ccs bam into fastq for minimap2 as the input file
$samtools fastq -@ 20 example.ccs.bam | gzip > example.ccs.fq.gz

# step 2, run minimap2, align the reads to the reference，here we used the hg19 which downloaded from UCSC.
minimap2 -Y -t 60 -a $mmi --cs --MD example.ccs.fq.gz | samtools view -Sb -@ 4 -o example.bam

# or align to the hg19 without indexed
# minimap2 -Y -t 60 -ax map-hifi $ref --cs --MD example.ccs.fq.gz | samtools view -Sb -@ 4 -o example.bam

# step 3, Split the bam file and generate a shell script to run
perl $sdir/ExtIDsGenShell.pl -i example.bam -d $dir/result -m example -n 2000 -r $ref -sn 1

# step4, Run the shell script, when search the overlap. We recommended to split the shell scrip and run them independently.
cd $dir/result/example
split -l 3 example.runFind.sh split.Find.
ls split.Find.* | grep -v log | perl -ne 'chomp;`nohup sh $_ &> $_.log &`;'
