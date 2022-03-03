#!/usr/bin/perl

# mycosnp_create_sample_sheet.pl v0.1.0
# Author: @hseabolt, @mjcipriano, @sateeshperi

# SYNOPSIS:
# Reads the contents of a given directory (expected to contain sequencing read data files in FASTQ format, optionally gzip-compressed),
# then uses standard file naming conventions to match up samples with their read mate information (including potentially multiple sequencing lanes),
# and finally, formats and outputs a csv sample sheet. 

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# GetOpts Variable declarations
my $input = "--";
my $output = "--";
my $version;
my $help;

sub usage {
	my $usage = "mycosnp_create_sample_sheet.pl v0.1.0\n
	PURPOSE: Reads the contents of a given directory (expected to contain sequencing read data files in FASTQ format, optionally gzip-compressed),
		then uses standard file naming conventions to match up samples with their read mate information (including potentially multiple sequencing lanes),
		and finally, formats and outputs a csv sample sheet. 

	\n
	USAGE:	mycosnp_create_sample_sheet.pl -i <input directory> -o <output>
	-i		input directory
	-o 		output file name
	-v	 	print version number and exit
	-h 		print this help message and exit
\n";
	print $usage;
}

GetOptions(	'input|i=s' => \$input, 
			'out|o=s'   => \$output,
			'version|v' => \$version,
			'help|h'    => \$help,
) or die usage();

# Print the version number or the help message and exit if -v or -h is activated
if ( $version ) 	{ die "mycosnp_create_sample_sheet.pl version 0.1.0\n"; }
if ( $help    )     { die usage();											}

################################################################################
# Read the input directory, capturing all files that end in .fastq (optionally .gz) or some flavor of this extension
opendir( my $dh, $input ) or die "mycosnp_create_sample_sheet::ERROR --> Cannot read the given input directory.  $!\n";
my @fastqs;
while ( readdir $dh )	{
	chomp $_;
	push @fastqs, $_ if ( $_ =~ /\.fq/ || $_ =~ /\.fastq/ );	# A curse upon ye who put .fastq somewhere in the filename other than the extension	
}

# For each FASTQ file:
# 1) deconstruct the string to get the sample/file name, lane info, and read information.
# 2) Then, populate a lookup table to store matching information per sample.
my %Samples = ();
foreach my $file ( @fastqs )	{
	
	# Deconstruct the sample name -- here I'm doing it with regex capture groups
	my ( $name, $fr_read ) = $file =~ /(.*)_[Rr]?(1|2)/;		
	my ( $lane ) = $name =~ s/_L(\d{3})?//;
	$lane = ( $1 eq $name )? "001" : $1;			# Default to 001 if there is no lane information
	
	# Add the file to the lookup table
	# Hierarchy is SAMPLENAME ==> LANE ==> READ (1|2) = FILENAME
	$Samples{$name}->{$lane}->{$fr_read} = $file;
}
	
# Set output filehandles
my $fileout = open( OUT, ">", "$output" ) if $output ne "--";
my $fhout;
if ( $fileout )		{	$fhout = *OUT;			}
else				{	$fhout = *STDOUT;		}	
	
# Print out a csv sample sheet with format <name>, <r1>, <r2>, <r1b>, <r2b>, ... <r1n>, <r2n>
foreach my $sample ( keys %Samples )		{
	my @line_out = ( $sample );
	foreach my $lane ( sort {$a <=> $b} keys %{$Samples{$sample}} )	{
		foreach my $read ( sort {$a <=> $b} keys %{$Samples{$sample}->{$lane}} )	{
			push @line_out, $Samples{$sample}->{$lane}->{$read};
		}
	}
	print $fhout join(",", @line_out), "\n";
}
close $fhout if ( $fileout );

exit;