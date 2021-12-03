#!/bin/bash

echo -n > array_design.txt

for i in {01..05} 
do 
    for j in {1..22}
    do 
	echo -e 04${i}'\t'${j} >> array_design.txt
    done
done
