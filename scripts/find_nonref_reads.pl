#!/usr/local/bin/perl -w
#
use strict;

my $Usage = qq!find_nonref_reads.pl <target bed file> <BAM file>\n!;

$#ARGV==1
    or die "$Usage";

my $target_bedfile = $ARGV[0];
my $bamfile = $ARGV[1];

my %element_ref_reads = ();
my %element_nonref_reads = ();
open REGIONS, $target_bedfile
    or die "Couldn\'t open $target_bedfile: $!\n";

while (<REGIONS>) {
    chomp;
    my ($chrom, $start, $end, $type, $refid, $rest) = split /\t/, $_;

    #print "$chrom:$start-$end:\n";

    open SAMREADS, "samtools view -F0x100 $bamfile $chrom:$start-$end | "
        or die "Couldn\'t run samtools: $!\n";

    while (<SAMREADS>) {
        chomp;
        my ($readname, $flag, $alignchrom, $alignstart, $score, $cigar, $restsam) = split /\t/, $_;

        my $left_clip = ($cigar =~ /^(\d+)[HS]/) ? $1 : 0;
        my $right_clip = ($cigar =~ /(\d+)[HS]$/) ? $1 : 0;

        my @large_indels = ();
        my @all_indels = ($cigar =~ /(\d+[DI])/g);
        my $no_indels = @all_indels;
        foreach my $thisindel (@all_indels) {
            my $size = ($thisindel =~ /(\d+)[DI]/) ? $1 : 0;
            if ($size >= 100) {
                push @large_indels, $size;
            }
        }
        if (@large_indels) {
            push @{$element_nonref_reads{$refid}}, $readname;
        }
        else {
            push @{$element_ref_reads{$refid}}, $readname;
        }
    }
    close SAMREADS; 

    my $no_ref = ($element_ref_reads{$refid}) ? @{$element_ref_reads{$refid}} : 0;
    my $no_nonref = ($element_nonref_reads{$refid}) ? @{$element_nonref_reads{$refid}} : 0;
    print "$refid\t$chrom:$start-$end\t$no_ref\t$no_nonref\n";
}
close REGIONS;

