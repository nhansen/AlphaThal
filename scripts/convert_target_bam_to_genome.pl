#!/usr/local/bin/perl -w
#

use File::Path qw( make_path );
use Getopt::Long;

use strict;

our $TOPDIR = $ENV{'LONGREADTOPDIR'};

if (!$TOPDIR) {
    die "You must set the LONGREADTOPDIR environment variable before running this script!\n";
}

our $SCRIPTSDIR=$TOPDIR.'/scripts';
#our $REFCOORDFILE='/data/nhansen/HERV_K_catalog/SVA_discovery/prep/sva_baits/ref_and_nonref_sva_regions.labeled.widerregion.bed';

my $Usage = qq!convert_target_bam_to_genome.pl <target bam file path> <target bed file>\n!;

$#ARGV==1
    or die "$Usage";

my $sample_bam = $ARGV[0];
my $target_bedfile = $ARGV[1];

my $rh_target_coords = read_target_coords($target_bedfile);

open SAM, "samtools view $sample_bam | ";

while (<SAM>) {
    chomp;
    my @fields = split /\t/, $_;

    my $target_id = $fields[2];
    my $pos = $fields[3];

    my $chrom = $rh_target_coords->{$target_id}->{chrom};
    my $target_start = $rh_target_coords->{$target_id}->{start};

    $fields[2] = $chrom;
    $fields[3] = $pos + $target_start;

    my $read_string = join "\t", @fields;
    print "$read_string\n";
}

close SAM;

sub read_target_coords {
    my $bedfile = shift;

    open REF, $bedfile
        or die "Couldn\'t open $bedfile: $!\n";

    my %target_info = ();
    while (<REF>) {
        chomp;
        my ($chrom, $start, $end, $target_id) = split /\t/, $_;
        $target_info{$target_id} = {chrom => $chrom,
                                  start => $start,
                                  end => $end};
    }
    close REF;

    return {%target_info};
}
