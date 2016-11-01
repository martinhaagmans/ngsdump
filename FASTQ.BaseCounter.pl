#!/bin/perl
$usage="\###
usage $0 extension

Deployed in a folder with read files with given extension gives %A %C %G %C per fastq file 
=======

###
" ;




$extension = $ARGV[0] ; 

if ($extension)	{	chomp $extension ; @files = <*$extension> ;

if (scalar(@files)==0)	{ @files = <*.$extension> ;  }

		}

@files = <*R1*.gz> unless ($extension) ; 


die $usage unless scalar(@files) != 0 ;


open(FH,">BasePercentages.txt") ;
print FH "Sample\tA\tT\tG\tC\n"  ;

foreach(@files) {	chomp $_ ;

@DNA= `zgrep -A1 "@" $_ | grep -v \"@\" | grep -v \"-\"` ;
$a = $b = $c = $d = 0 ;

foreach $DNA (@DNA) 	{

$a+=($DNA=~tr/A//);
$b+=($DNA=~tr/C//);
$c+=($DNA=~tr/G//);
$d+=($DNA=~tr/T//);	
			}
$Total=$a+$b+$c+$d;	

$aperc = $a/$Total ;
$bperc = $b/$Total ;
$cperc = $c/$Total ;
$dperc = $d/$Total ;

print FH "$_\t$aperc\t$dperc\t$cperc\t$bperc\n"  ;
print "$_\t$aperc\t$dperc\t$cperc\t$bperc\n"  ;

		}

close FH ;


exit;
