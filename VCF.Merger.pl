#!/usr/bin/perl

use strict ;
use warnings ;
use Config::Simple  ;
use Cwd 'abs_path';

my $usage="\###
usage $0

Deployed in a folder with VCF files with combines all VCF files into combined.vcf
=======

###
" ;

############read config file

my $config = "/home/corona/Desktop/martin/ZaaGnome.INI" ; my %config ; 
Config::Simple->import_from ("$config", \%config ) ;
my($reference,$java,$gatk) ;  
my @files ; 

$java=$config{"SYSTEM.java7" } ; 
$reference=$config{"ALIGNMENT.reference" } ; 
$gatk=$config{"GATK.tp53gatk" } ; 
my $gatk4 = "$java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar" ;

my @vcf = <*.vcf> ;

open (FH,">script.sh") ;
print FH "\#!/bin/sh\n" ;
print FH "$gatk4 -T CombineVariants  -R $reference    -filteredRecordsMergeType   KEEP_IF_ANY_UNFILTERED \\\n" ;

foreach (@vcf)	{ chomp $_ ; print FH "-V $_ \\\n"	; } 

print FH " -o combined.vcf \n"; 

execute("bash script.sh && rm scipt.sh") ;
exit;

##############################################################################

sub execute {
 my $syscall = shift;
 my $sysreturn = system($syscall);
 die "$syscall WENT WRONG\n" if $sysreturn > 0 ; 
} 	
