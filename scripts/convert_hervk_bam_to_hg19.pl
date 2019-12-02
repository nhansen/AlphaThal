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
our $HEADERFILE='/fdb/genome/hg19/hg19.fa.fai';
our $REFCOORDFILE='/data/nhansen/HERV_K_catalog/SVA_discovery/prep/sva_baits/ref_and_nonref_sva_regions.labeled.widerregion.bed';

my $Usage = qq!convert_hervk_bam_to_hg19.pl <HERVK bam file path>\n!;

$#ARGV==0
    or die "$Usage";

my $sample_bam = $ARGV[0];

open HEADER, $HEADERFILE
    or die "Couldn\'t open header file $HEADERFILE: $!\n";

while (<HEADER>) {
    #print;
}
close HEADER;

my $rh_hervk_coords = read_hervk_coords();

open SAM, "samtools view $sample_bam | ";

while (<SAM>) {
    chomp;
    my @fields = split /\t/, $_;

    my $hervk_id = $fields[2];
    my $pos = $fields[3];

    my $chrom = $rh_hervk_coords->{$hervk_id}->{chrom};
    my $hervk_start = $rh_hervk_coords->{$hervk_id}->{start};

    $fields[2] = $chrom;
    $fields[3] = $pos + $hervk_start;

    my $read_string = join "\t", @fields;
    print "$read_string\n";
}

close SAM;

sub read_hervk_coords {

    open REF, $REFCOORDFILE
        or die "Couldn\'t open $REFCOORDFILE: $!\n";

    my %hervk_info = ();
    while (<REF>) {
        chomp;
        my ($chrom, $start, $end, $hervk_id) = split /\t/, $_;
        $hervk_info{$hervk_id} = {chrom => $chrom,
                                  start => $start,
                                  end => $end};
    }
    close REF;

    #open NONREF, $NONREFCOORDFILE
        #or die "Couldn\'t open $NONREFCOORDFILE: $!\n";
#
    #while (<NONREF>) {
        #chomp;
        #my ($chrom, $start, $end, $hervk_id) = split /\t/, $_;
        #$hervk_info{$hervk_id} = {chrom => $chrom,
                                  #start => $start,
                                  #end => $end};
    #}
    #close NONREF;
    return {%hervk_info};
}
