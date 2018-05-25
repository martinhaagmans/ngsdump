#!/bin/bash

DIR=$(dirname $0)

for i in $(ls $DIR/MS*.txt) ; do
	pyCNV --create -c BLA -i $i && mv $i $DIR/toegevoegd
done ;


