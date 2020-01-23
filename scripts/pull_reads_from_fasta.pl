#!/usr/local/bin/perl -w
#
use strict;

use FileHandle;

my $Usage = "pull_reads_from_fasta.pl <file of readnames> <fasta file>\n";

$#ARGV==1
    or die "$Usage";

my $readfile = $ARGV[0];
my $fastafile = $ARGV[1];

my %readnames = ();

my $readfilestring = ($readfile =~ /\.gz$/) ? "gunzip -c $readfile | " : $readfile;
my $read_fh = FileHandle->new($readfilestring);

while (<$read_fh>) {
    if (/^\s*(\S+)\s*/) {
        $readnames{$1} = 1;
    }
}
close $read_fh;

my $fastafilestring = ($fastafile =~ /\.gz$/) ? "gunzip -c $fastafile | " : $fastafile;

my $fa_fh = FileHandle->new($fastafilestring);

my $printthisread = 0;
while (<$fa_fh>) {
    my $thisline = $_;
    if ($thisline =~ /^>\s*(\S+)/) { # correctly formatted entry line
        if ($readnames{$1}) { # interested in this read
            $printthisread = 1;
            delete $readnames{$1};
            print $thisline;
        }
        else {
            $printthisread = 0;
        }
    }
    else {
        if ($printthisread) {
            print $thisline;
        }
    }
}
close $fa_fh;

# print remaining reads to stderr:
#
foreach my $forgottenread (keys %readnames) {
    print STDERR "UNFOUND\t$forgottenread\n";
}
