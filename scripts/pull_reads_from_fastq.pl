#!/usr/local/bin/perl -w
#
use strict;

use GTB::File qw(Open);

my $Usage = "pull_reads_from_fastq.pl <file of readnames> <fastq file>\n";

$#ARGV==1
    or die "$Usage";

my $readfile = $ARGV[0];
my $fastqfile = $ARGV[1];

my %readnames = ();
my $read_fh = Open($readfile);

while (<$read_fh>) {
    if (/^\s*(\S+)\s*/) {
        $readnames{$1} = 1;
    }
}
close $read_fh;

my $fq_fh = Open($fastqfile);

my $printthisread = 0;
while (<$fq_fh>) {
    if (/^@\s*(\S+)/) { # correctly formatted first line
        if ($readnames{$1}) { # interested in this read
            $printthisread = 1;
            delete $readnames{$1};
        }
        else {
            $printthisread = 0;
        }
        if ($printthisread) {
            print;
        }
        for (my $nextline = 1; $nextline <= 3; $nextline++) {
            if (my $line = <$fq_fh>) {
                if ($printthisread) {
                    print $line;
                }
            }
            else {
                die "Unexpected formatting in fastq file: must have four lines per entry!\n";
            }
        }
    }
}
close $fq_fh;

# print remaining reads to stderr:
#
foreach my $forgottenread (keys %readnames) {
    print STDERR "UNFOUND\t$forgottenread\n";
}
