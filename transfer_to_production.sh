#!/bin/bash

DIR=$(pwd) ;


for i in $(ls $DIR/MS*.txt) ; do
	scp $i ux-p-dnadiag:$DIR ;
	rm -f $i ;
done ;

ssh ux-p-dnadiag "bash $DIR/add_to_db.sh" 

