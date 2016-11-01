#!/usr/bin/perl -w

use strict ;

my $usage = "$0 Genes.txt CaptureTarget.bed \n"  ;

my $out = 'CoreTest.bed' ;

die $usage unless scalar(@ARGV) == 2 ;

`rm $out && touch $out` ;

open (FH,$ARGV[0]) ; my @genes = <FH> ; close FH ;



foreach (@genes)    {

    chomp $_ ;

     `awk \'(\$4==\"$_\") {print}\' $ARGV[1]  >> $out  `  ;

     my $nrtargets = `awk \'(\$4==\"$_\") {print}\' $ARGV[1] | wc -l`  ;

     chomp $nrtargets ;

     if ($nrtargets == 0)   { print $_ . "\n"  ;}


}


exit ;
