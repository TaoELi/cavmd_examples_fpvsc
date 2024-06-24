#!/bin/bash

filename=init.xyz

touch $filename

for i in {1..36}
do
    cat single_cell.xyz >> $filename
done

