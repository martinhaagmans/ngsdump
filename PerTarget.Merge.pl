#!/usr/bin/perl

use strict ;
use warnings ;

my $usage="\###
usage $0 targetBED extension

Deploy in a folder with Picard's PerTargetCoverage files with targetBED \
Combines all Relative Coveage into file RelativeMeanCoverage.txt and all \
Mean Coverage into file MeanCoverage.txt
=======

###
" ;

my @files ;
my $capture = $ARGV[0] ;  die $usage  unless $capture ;

my $extension = $ARGV[1] ;

if ($extension)	{
    chomp $extension ; @files = <*$extension> ;
    if (scalar(@files)==0)	{
        @files = <*.$extension> ;
    }
}

@files = <*.PerTarget>  unless ($extension) ;

die $usage unless scalar(@files) != 0 ;

open (FHOUT, ">MeanCoverage.txt") ;
open (FHRELOUT, ">RelativeMeanCoverage.txt") ;

open (FH, "<$capture") ; my @regions = <FH> ; close FH ;

foreach my $region (@regions) {
    chomp $region ;
    $region =~ s/^\s+//;
    $region =~ s/\s+$//;
    print FHOUT "$region\t" ;
    print FHRELOUT "$region\t" ;

    my ($targetchrom, $targetstart, $targetend) = split("\t", $region) ;


    if ($targetchrom eq "chrom") {
        foreach (@files)	{
            chomp $_ ;
            print FHOUT  "$_\t\t" ;
            print FHRELOUT  "$_\t\t" ;
        }

        print FHOUT "\n" ;
        print FHRELOUT "\n" ;
        next ;
    }

    foreach (@files)	{
        chomp $_ ;
        open (FH, "<$_") ; my @coverage = <FH> ; close FH ;

        foreach my $line (@coverage) {
            chomp $line ;
            $line =~ s/\s+/\t/;
            my ($chrom, $start, $end, $size, $capture, $gcperc, $meancov, $normalizedcov) = split("\t", $line) ;
            $meancov =~ s/\./,/g ;
            $normalizedcov =~ s/\./,/g ;

            if ($targetchrom eq $chrom && $targetstart == $start && $targetend == $end) {
                print FHOUT  "$meancov\t" ;
                print FHRELOUT  "$normalizedcov\t" ;
            }

        }
    }

    print FHOUT "\n" ;	print FHRELOUT  "\n" ;
}

close FHOUT ; close FHRELOUT ;

    exit ;
