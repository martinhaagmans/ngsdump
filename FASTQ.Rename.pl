#!/usr/bin/perl
#hernoemt fastq bestanden aan de hand van een lijst met 2 kolommen, eerste kolom wordt vervangen door tweede kolom
$input = $ARGV[0] ;
open (FH,"<$input") ; @input  = <FH> ; close FH ;

foreach (@input) 	{ chomp $_ ; ($nr,$dnr) = split (/\t/,$_) ; `rename $nr\_S $dnr\_S *` ; } 

exit ;



sub execute {
 my $syscall = shift;
 my $sysreturn = system($syscall);
 die "$syscall WENT WRONG\n" if $sysreturn > 0 ; 
} 	
