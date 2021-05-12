#!/bin/bash
#SBATCH --cpus-per-task 16
fidx=1
while(( $fidx<=70 ))
do
	echo $fidx
	nohup matlab -nosplash -nodesktop -r "fidx=$fidx;ln_xcorr_bz;quit"
	let "fidx++"
done
