#!/bin/bash -l
#insert use commands here
use Java-1.8

#  $1 is the task file

# input the $SGE_TASK_IDth line of this file into $id
export id=`awk "NR==$SGE_TASK_ID" $1`
#split id on whitespace into array that can be accessed like ${splitID[0]}
#splitID=($id)
#export id=${splitID[0]}
echo  $SGE_TASK_ID : $id
#redirect stdout and stderr to file
export logDir="logs"
mkdir -p $logDir
exec 1>>$logDir/$id.olog
exec 2>&1
set -e

#output directory for alignments
export alignDir="aligned"
mkdir -p $alignDir

if [ ! -e $logDir/$id.done ]
then
	echo Job $JOB_ID:$SGE_TASK_ID started on $HOST: `date`
	outid=$alignDir/${id/paired.fastq/align.txt}
	corrid=${outid/align.txt/corr_align.txt}
	java -Xmx16G -jar /broad/ragon_compute/scripts/TCR_BCR_alignment/migmap-1.0.2/migmap-1.0.2.jar --blast-dir /broad/ragon_compute/scripts/igblast-bin/binaries/linux-x64 -R IGH -S human $id $outid
	java -Xmx16G -cp /broad/ragon_compute/scripts/TCR_BCR_alignment/migmap-0.9.9/migmap-1.0.0-SNAPSHOT.jar com.antigenomics.migmap.Correct $outid $corrid

	touch $logDir/$id.done
	echo Job finished: `date`
else
	echo Job already done. To redo, rm $logDir/$id.done
fi
