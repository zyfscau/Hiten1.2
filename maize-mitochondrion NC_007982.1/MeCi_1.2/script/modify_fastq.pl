#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Command line options
my $input_file;
my $output_file;
my $reverse = 0; # Boolean, default is false (no reverse)
my $read_number = 1; # Default read number

# Get options
GetOptions(
    "input|i=s" => \$input_file,
    "output|o=s" => \$output_file,
    "reverse|r!" => \$reverse, # Boolean flag
    "read|rn=i" => \$read_number
) or die "Error in command line arguments\n";

# Check mandatory arguments
unless ($input_file && $output_file) {
    die "Usage: perl script.pl --input [file] --output [file] [--reverse] [--read 1|2]\n";
}

# Check read number validity
unless ($read_number == 1 || $read_number == 2) {
    die "Read number must be 1 or 2\n";
}

# Open the input and output FASTQ files
open(my $in, "<", $input_file) or die "Cannot open $input_file: $!";
open(my $out, ">", $output_file) or die "Cannot open $output_file: $!";

while (my $header = <$in>) {
    # Process the header line
    my ($read_name, $rest) = split(/\s/, $header, 2);
    $read_name .= "_$read_number";
    print $out "$read_name\n";

    # Process the sequence line
    my $seq = <$in>;
    chomp($seq);
    if ($reverse) {
        $seq = reverse $seq;
        $seq =~ tr/ACGTacgt/TGCAtgca/;
    }
    print $out "$seq\n";

    # Process the '+' line
    my $plus = <$in>;
    print $out $plus;

    # Process the quality line
    my $qual = <$in>;
    chomp($qual);
    if ($reverse) {
        $qual = reverse $qual;
    }
    print $out "$qual\n";
}

close($in);
close($out);
