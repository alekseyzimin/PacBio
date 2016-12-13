#!/bin/bash
#this script runs delta-filter in parallel, input is DELTAFILE.delta, output to STDOUT
DELTAFILE=$1
OPTIONS=$2
NUM_THREADS=$3
PID=$$;

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

set -e
#run delta-filter commands
COMMAND="head -n 2 $DELTAFILE.delta > $PID.head && tail -n +3 $DELTAFILE.delta | ufasta split "
for i in $(seq 1 $NUM_THREADS);do
    COMMAND=${COMMAND}" >(cat $PID.head /dev/stdin |delta-filter $OPTIONS /dev/stdin | tail -n +3 > $PID.$i.fdelta && touch $PID.$i.success || touch $PID.$i.fail) "
done
eval $COMMAND

#wait for subshells
DONE=0
COUNTER=0
until [ $COUNTER -eq  $NUM_THREADS ];do
    sleep 1
    COUNTER=0
    for i in $(seq 1 $NUM_THREADS);do
	if [ -e $PID.$i.success ];then 
	    let COUNTER+=1
	fi
	if [ -e $PID.$i.fail ];then
	    echo "delta-filter job $i failed, exiting"
	    for i in $(seq 1 $NUM_THREADS);do
		rm -f $PID.$i.fdelta $PID.$i.success $PID.$i.fail
	    done
	    exit 1
	fi
    done
done

#concatenate outputs
COMMAND="cat $PID.head "
for i in $(seq 1 $NUM_THREADS);do
    COMMAND=${COMMAND}"  $PID.$i.fdelta "
done
COMMAND=${COMMAND}" | delta-filter $OPTIONS /dev/stdin > $DELTAFILE.fdelta"
eval $COMMAND

#cleaning up
rm -f $PID.head
for i in $(seq 1 $NUM_THREADS);do
    rm -f $PID.$i.fdelta $PID.$i.success
done
