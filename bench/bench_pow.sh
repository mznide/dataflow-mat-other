#!/bin/bash


EXE=$1

name=${EXE#../}
name=${name%%/*}
logfile="$name.txt"


test -f $logfile && exit

m=0
while (( m <= 2)); do
n=3072
p=1
	while (( p <= 30 )); do
		sleep 2
		$EXE --size $n -p $p | tee -a $logfile
		(( p +=1 ))
	done

	(( m += 1 ))
done
