#!/bin/bash

for i in {1..15}
do
	eps=$(echo "0.1*$i" | bc)
	echo $eps | Executables/initialConditions
	echo $eps | Executables/evolution
done



