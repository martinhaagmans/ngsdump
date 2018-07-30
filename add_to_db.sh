#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

for i in $(ls $DIR/MS*.txt) ; do
	CNV --create -c BLA -i $i && mv $i $DIR/toegevoegd
done ;


