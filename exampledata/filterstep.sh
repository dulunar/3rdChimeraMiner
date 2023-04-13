#!/bin/bash
#########################################################################
# FileName: filterstep.sh
# Version: c0c175c6-1dee-4559-b041-933ec4c64f7d
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Thu 13 Apr 2023 12:34:53 AM CST
#########################################################################

dir=`pwd`
sdir=`dirname $dir`
samp=example

while :
do
	psi=`ps xu | grep perl | grep FindChimeras.latest.pl | grep ${samp} | tail -1`
	if [[ $psi ]]
	then
		sleep 600
	else
		# Downstream analysis and achieve valid chimeras
		# Here, we assume that the length threshold of the overlapping sequence was set to 2 nt.
		perl $sdir/ChimerasDownstream.pl -d $dir/result/$samp -od $dir/result/$samp/chimeras_analysis_lap2 -n $samp -p 2 

		# Split chimeric read into multiple normal segments
		perl $sdir/SplitChimericRead.pl -i $dir/$samp.ccs.bam -f $dir/result/$samp/chimeras_analysis_lap2/$samp.chimera.id -c $dir/result/$samp/chimeras_analysis_lap2/$samp.normal.txt -o $dir/result/$samp/chimeras_analysis_lap2/$samp.chimera.fq

		break
	fi
done
