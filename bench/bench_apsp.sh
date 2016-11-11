#!/bin/bash


EXE=$1

name=${EXE#../}
name=${name%%/*}
logfile="$name.txt"


test -f $logfile && exit

m=0
while (( m <= 1)); do

n=512
	while (( n <= 3072 )); do
		sleep 2
		$EXE --size $n | tee -a $logfile
		(( n +=512 ))
	done
	(( m += 1))
done
