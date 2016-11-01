#!/usr/bin/perl

use strict ;
use warnings ;
use Config::Simple  ;


my $usage="\###
usage $0 VCF ID-list

Outputs Select.out.vcf with variants from VCF from all ID's in ID-list
=======

###
" ;

############read config file

my $scriptdir = abs_path($0) ; $scriptdir =~ s/\/\w*\.pl//gi ;  my $config = "$scriptdir../ZaaGnome.INI" ; my %config ;
Config::Simple->import_from ("$config", \%config ) ;
my($reference,$java,$gatk) ;  


$java=$config{"SYSTEM.java7" } ; 
$reference=$config{"ALIGNMENT.reference" } ; 
$gatk=$config{"GATK.wesgatk" } ; 

my $gatk4 = "$java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar" ;

my $inputfile = $ARGV[0] ;
my $list = $ARGV[1] ;

die $usage unless ($inputfile && $list) ;

open (FH,"<$list") ; my  @vcf = <FH> ; close FH ;

open (FH,">script.sh") ;
print FH "\#!/bin/sh\n" ;
print FH "$gatk4 -T SelectVariants  -R $reference  -o Select.out.vcf \\\n" ;

foreach (@vcf)	{ chomp $_ ; print FH "-sn $_ \\\n"	; } 

print FH "-V $inputfile -env -nt 6 \n"; 

execute("bash script.sh && rm script.sh") ;

exit;

##############################################################################

sub execute {
 my $syscall = shift;
 my $sysreturn = system($syscall);
 die "$syscall WENT WRONG\n" if $sysreturn > 0 ; 
} 	
