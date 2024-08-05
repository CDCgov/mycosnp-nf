#!/usr/bin/env perl

# mycosnp_combine_lanes.pl v0.1.0
# Author: @hseabolt, @mjcipriano, @sateeshperi

# SYNOPSIS:
# Reads and parses a CSV sample sheet from mycosnp_create_sample_sheet.pl to combine lanes with system calls to 'cat',
# Then creates a new CSV sample sheet with the fields <name>, <r1>, <r2>

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# GetOpts Variable declarations
my $input = "--";
my $output = "--";
my $version;
my $help;

sub usage {
	my $usage = "\nmycosnp_combine_lanes.pl\nv0.1.0\n
	PURPOSE: Reads and parses a CSV sample sheet from mycosnp_create_sample_sheet.pl to combine lanes with system calls to 'cat',
		 Then creates a new CSV sample sheet with the fields <name>, <r1>, <r2>

	USAGE:	mycosnp_combine_lanes.pl -i <input sample sheet> -o <output sample sheet>
	
	ARGUMENTS:
	-i | --input		DIR (Required).  Input sample sheet, expecting csv format.
	-o | --output		STR (Optional).  Output sample sheet file name.  If no argument given, prints to STDOUT.	
	-v | --version	 	Print version number and exit.
	-h | --help		Print this help message and exit.
	
	\n";
	print $usage;
}

GetOptions(	'input|i=s'  => \$input, 
			'out|o=s'    => \$output,
			'version|v'  => \$version,
			'help|h'     => \$help,
) or die usage();

# Print the version number or the help message and exit if -v or -h is activated
if ( $version ) 	{ die "mycosnp_combine_lanes.pl version 0.1.0\n"; }
if ( $help    )     { die usage();											}

################################################################################
# Read the sample sheet, expecting a csv format
my $fh = *STDIN;
my $succin = open(CSV, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *CSV if ( $succin ); 

# Set output filehandles
my $succout = open( OUT, ">", "$output" ) if $output ne "--";
my $fhout;
if ( $succout )		{	$fhout = *OUT;			}
else				{	$fhout = *STDOUT;		}

# Print sample sheet header line
print $fhout join(",", "sample", "fastq_1", "fastq_2" ), "\n";

while ( <$fh> )		{
	chomp $_;
	my @line = split(",", $_);
	my $seqid = shift @line;
	my $ext = "fastq.gz"; 			# Defaulting to .gz extension
	
	system("mkdir $seqid");
	
	if ( scalar @line > 3 )	{
		$ext = ( $line[1] =~ /\.gz$/ )? "fastq.gz" : "fastq";
		for ( my $i=1; $i < scalar @line; $i++ )	{
			my $reads = ( $i % 2 == 0 )? "R2" : "R1";	
			system("cat $line[$i] >> $seqid/$seqid\_$reads.$ext");
		}
	}
	else 	{
		$ext = ( $line[2] =~ /\.gz$/ )? "fastq.gz" : "fastq";
		system("ln -s $line[2] $seqid/$seqid\_R1.$ext");
		system("ln -s $line[2] $seqid/$seqid\_R2.$ext");
	}
	
	# Print entry out to the sample sheet with the qualified full path
	my $realpathR1 = `realpath $seqid/$seqid\_R1.$ext`; chomp $realpathR1;
	my $realpathR2 = `realpath $seqid/$seqid\_R2.$ext`; chomp $realpathR2;
	print $fhout join(",", $seqid, $realpathR1, $realpathR2), "\n";
}
close $fh if ( $succin );
close $fhout if ( $succout );

exit;