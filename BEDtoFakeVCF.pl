#!/usr/bin/perl

use strict ;
use warnings ;

my $usage="\###
usage $0 BED-file
creates fake vcf for every position in BED
=======

###
" ;

my $bed = $ARGV[0] ; die $usage unless $bed ; open(FH,"<$bed") ;

my $header = "\##fileformat=VCFv4.1
##FILTER=<ID=LowQual,Description=\"Low quality\">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##contig=<ID=chrM,length=16571,assembly=hg19>
##reference=file:///home/corona/gatkrecources/lifescope.hg19.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CrashDummy" ;

my $restvcf = ".	C	G	2227.77	PASS	AC=2;AF=1.00;AN=2;DP=71;	GT:AD:GQ:PL	1/1:0,71:99:2256,213,0" ;


##All bed positions
open(OUT,">$bed.vcf") ;
print OUT "$header\n" ;
while (<FH>)	{ chomp ;
my($chr,$start,$end) = split /\t/ ;

	while ($start <=$end)	{
	print OUT "$chr\t$start\t$restvcf\n" ; $start++ ;

	}

next ;

}


close FH ;
close OUT ;


exit ;
