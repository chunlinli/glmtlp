#!/bin/bash
for round in 1 2 3
do
    for i in 4 3 2 1
    do
	echo "Running round $round test $i..."
	R CMD BATCH test$i.R
    done
done
