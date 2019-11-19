#!/usr/local/bin/perl -w
#

use Getopt::Long;

use strict;

our $SCRIPTSDIR='/data/nhansen/HERV_K_catalog/discovery2/scripts';
our $REFCOORDFILE='/data/nhansen/HERV_K_catalog/discovery2/prep/allele_sequences/herv_k_hml2_rm_coords.merged.5000.ge600.widerregion.bed';
our $NONREFCOORDFILE='/data/nhansen/HERV_K_catalog/discovery2/prep/allele_sequences/herv_k_hml2_nonref_coords.widerregion.bed';
our $ELEMENTDIR='/data/nhansen/HERV_K_catalog/discovery2/prep/hervk_elements';

my $Usage = qq!create_annotation_bedfile.pl\n!;

$#ARGV<0
    or die "$Usage";

my $rh_hervk_coords = read_hervk_coords();

opendir ELEMENTS, $ELEMENTDIR
    or die "Couldn\'t open $ELEMENTDIR for reading: $!\n";

my @cmfiles = grep /\.cm\.out$/, readdir ELEMENTS;

closedir ELEMENTS;

foreach my $cm_file (@cmfiles) {
    open CM, "$ELEMENTDIR/$cm_file"
        or die "Couldn\'t open $ELEMENTDIR/$cm_file: $!\n";

    while (<CM>) {
        if (/^ALIGNMENT\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/) {
            my ($score, $hervk_id, $element, $start, $end) = ($1, $5, $9, $6, $7);
            my $chrom = $rh_hervk_coords->{$hervk_id}->{'chrom'};
            my $hervk_start = $rh_hervk_coords->{$hervk_id}->{'start'};

            my $hervk_end = $hervk_start + $end - 1;
            $hervk_start += $start - 1;

            print "$chrom\t$hervk_start\t$hervk_end\t$element\t$score\t+\n";
        }
        elsif (/^ALIGNMENT\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+C\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/) {
            my ($score, $hervk_id, $element, $start, $end) = ($1, $5, $9, $6, $7);
            my $chrom = $rh_hervk_coords->{$hervk_id}->{'chrom'};
            my $hervk_start = $rh_hervk_coords->{$hervk_id}->{'start'};

            my $hervk_end = $hervk_start + $end - 1;
            $hervk_start += $start - 1;

            print "$chrom\t$hervk_start\t$hervk_end\t$element\t$score\t\-\n";
        }
    }

    close CM;
}


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

    open NONREF, $NONREFCOORDFILE
        or die "Couldn\'t open $NONREFCOORDFILE: $!\n";

    while (<NONREF>) {
        chomp;
        my ($chrom, $start, $end, $hervk_id) = split /\t/, $_;
        $hervk_info{$hervk_id} = {chrom => $chrom,
                                  start => $start,
                                  end => $end};
    }
    close NONREF;
    return {%hervk_info};
}
