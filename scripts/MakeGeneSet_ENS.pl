#!/usr/bin/perl
use strict;
use warnings;


if (@ARGV <2)
{
    print "
    This script takes an extracted gene info TF_to_gene file and generates a gmt geneset file for use with GSEA, using Ensembl IDs.

    REQUIRED:
    -i input extracted gene info file
    -o output gmt file\n";
	exit;
}

my $input="NA";
my $output=".";

my $n=-1;
foreach(@ARGV)
{
	$n++;
	if($_=~/-i/) {$input=$ARGV[$n+1]; next;}
	if($_=~/-o/) {$output=$ARGV[$n+1]; next;}
}
if (! -e $input){print "ERROR: Your input file ($input) does not exist\n";exit;}


open (FH, "<$input");
open (OUT, ">$output");
my $TF="NA";
my @array;
while(defined(my $line=<FH>))
{
	chomp $line;
	my @info=split(/\t/,$line);
	my $ens=$info[1];
	chomp $ens;
	my @ENS=split(/\./,$ens);
	if ($info[0] eq $TF)
	{
		#print OUT "\t$ENS[0]";
		push (@array, $ENS[0]);
	}
	else
	{
		if ($TF ne "NA")
		{
			print OUT "$TF\tNA\t";
			my @filtered = uniq(@array);
			print OUT  join( "\t", @filtered ), "\n";
			@array=();
		}
		$TF=$info[0];
		#print OUT "$TF\tNA\t$ENS[0]";
		push (@array, $ENS[0]);
	}
}
close FH;
close OUT;

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
