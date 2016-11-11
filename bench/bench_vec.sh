#!/bin/bash


EXE=$1

name=${EXE#../}
name=${name%%/*}
logfile="$name.txt"


test -f $logfile && exit

m=0
while (( m <= 2)); do
n=3072
vec_nu=1024
	while (( vec_nu <= 20480 )); do
		sleep 2
		$EXE --size $n -v $vec_nu | tee -a $logfile
		(( vec_nu += 1024 ))
	done

	(( m += 1 ))
done
