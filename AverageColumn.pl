#!/usr/bin/perl

use strict ;
use warnings ;

my $usage="\###
usage $0 extension column

Deploy in a folder to get the average of a specific column of all files with extension
###
" ;

my $extension = $ARGV[0] ; die "$usage\n\nERROR: no extension\n\n" unless $extension ;
my $column    = $ARGV[1] ; die "$usage\n\nERROR: no column\n\n" unless $column ;

my @files = <*$extension> ;

if (scalar(@files)==0)	{	@files = <*.$extension> ;	}


die $usage unless scalar(@files) != 0 ;

open(FH," >Average.$column.out") ;
foreach (@files)	{ next if ($_ eq "Average.$column.out") ;


my @array = `awk \' { print \$$column} \' $_` ; chomp @array ; my $header = shift(@array) ;

my $ave = &average(\@array);
my $std = &stdev(\@array);

print "$_\t$ave\t$std\n" ;
print FH "$_\t$ave\t$std\n" ;

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

			}
