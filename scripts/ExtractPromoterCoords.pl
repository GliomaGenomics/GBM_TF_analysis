#!/usr/bin/perl 

use strict;
use warnings;

if (@ARGV <2)
{
    print "
    This script takes a gtf file, identifies the TSS fpr each transcript and grabs
    the coords 1000bp either side to attempt to capture the promoter. The output will be
    the same name as the gtf but with .gtf replaced with _PromotersTSS_size.txt

    REQUIRED:
    -g [gtf] {MUST use full filepath}
    -s [size] of promoter i.e. bp up and downstream to extract {default:1000}
    -o [outputDirectory] {default:pwd}\n";
	exit;
}

my $gtf="NA";
my $outDir=".";
my $size=1000;

my $n=-1;
foreach(@ARGV)
{
	$n++;
	if($_=~/-g/) {$gtf=$ARGV[$n+1]; next;}
	if($_=~/-o/) {$outDir=$ARGV[$n+1]; next;}
	if($_=~/-s/) {$size=$ARGV[$n+1]; next;}
}
if ($gtf eq "NA"){print "ERROR: You must supply a gtf file using -g\n";exit;}
if (! -e $gtf){print "ERROR: Your input file ($gtf) does not exist\n";exit;}

my @paths=split(/\//,$gtf);
my $resFile=pop(@paths);
my $nameString="_PromotersTSS_".$size;
$resFile=~s/\.gtf/$nameString\.txt/;
open (RES, ">$outDir/$resFile");
print RES "chrom\tstart\tend\tstrand\tGeneID\tGeneName\tTranscriptID\tTranscriptName\tType\n";
my $chr;
my $tss;
my $start;
my $end;
my $strand;
my $geneID;
my $geneName;
my $transcriptID;
my $transcriptName;
my $type;

open (FH,"<$gtf") or die "unable to open $gtf\n";
while(defined(my $line=<FH>))
{
	if ($line=~/^#/){next;}
	chomp $line;
	my @gtfdata=split(/\t/,$line);
	if ($gtfdata[2] ne "transcript"){next;}
	if ($gtfdata[6] eq "+"){$tss=$gtfdata[3];}
	else {$tss=$gtfdata[4];}
	if ($tss<1001){$start=1;}
	else{$start=$tss-$size;}
	$end=$tss+$size;
	$strand=$gtfdata[6];
	$chr=$gtfdata[0];
	if ($chr eq "chrM"){$chr="chrMT";} #This is needed because GTRD uses chrMT
	my @ids=split(/"/,$gtfdata[8]);
	$geneID=$ids[1];
	$geneName=$ids[7];
	$transcriptID=$ids[3];
	$transcriptName=$ids[11];
	$type=$ids[5];
	print RES "$chr\t$start\t$end\t$strand\t$geneID\t$geneName\t$transcriptID\t$transcriptName\t$type\n";
}
