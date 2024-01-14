#!/bin/bash
for i in $(seq -f "%03g" 1 25)
do
	echo $i
	mkdir traj${i}
done

