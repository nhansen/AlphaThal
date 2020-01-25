#!/usr/local/bin/perl -w
#
use strict;

my $Usage = qq!create_read_table.pl <target bed file> <minimap2 bam file> <MASHmap output file>\n!;

$#ARGV==2
    or die "$Usage";

my $target_bedfile = $ARGV[0];
my $minimap_bam = $ARGV[1];
my $mashmap_file = $ARGV[2];

my $rh_mash_hits = read_mash_hits($mashmap_file);

open REGIONS, $target_bedfile
    or die "Couldn\'t open $target_bedfile: $!\n";

print "RefID\tRead\tReadLength\tLeftMashLength\tLeftMashIdent\tLeftClip\tRightMashLength\tRightIdent\tRightClip\n";

while (<REGIONS>) {
    chomp;
    my ($chrom, $start, $end, $type, $refid, $rest) = split /\t/, $_;

    open SAMREADS, "samtools view -F0x100 $minimap_bam $chrom:$start-$end | "
        or die "Couldn\'t run samtools: $!\n";
   
    while (<SAMREADS>) {
        chomp;
        my ($readname, $flag, $alignchrom, $alignstart, $score, $cigar, $restsam) = split /\t/, $_;
  
        my $left_clip = ($cigar =~ /^(\d+)[HS]/) ? $1 : 0;
        my $right_clip = ($cigar =~ /(\d+)[HS]$/) ? $1 : 0;
       
        my $ra_mash_hits = ($rh_mash_hits->{$readname}) ? $rh_mash_hits->{$readname} : [];
        my @lefthits = grep {$_->{chrom} eq $alignchrom && $_->{regionend} == $start} @{$ra_mash_hits};
        my $mashlength_left = (@lefthits) ? $lefthits[0]->{hitlength} : 'NA';
        my $mashident_left = (@lefthits) ? $lefthits[0]->{perciden} : 'NA';
        my @righthits = grep {$_->{chrom} eq $alignchrom && $_->{regionstart} == $end} @{$ra_mash_hits};
        my $readlength = (@lefthits) ? $lefthits[0]->{readlength} : ((@righthits) ? $righthits[0]->{readlength} : 'NA');
        my $mashlength_right = (@righthits) ? $righthits[0]->{hitlength} : 'NA';
        my $mashident_right = (@righthits) ? $righthits[0]->{perciden} : 'NA';
        print "$refid\t$readname\t$readlength\t$mashlength_left\t$mashident_left\t$left_clip\t$mashlength_right\t$mashident_right\t$right_clip\n";
    }
    close SAMREADS; 
}
close REGIONS;

sub read_mash_hits {
    my $mash_file = shift;

    my %mash_hits = ();
    open MASH, $mash_file
        or die "Couldn\'t open $mash_file for reading: $!\n";
    while (<MASH>) {
        chomp;
        my ($readname, $readlength, $readstart, $readend, $readstrand, $refhit, $reflength, $refstart, $refend, $percident) = split /\s/, $_;

        my ($chrom, $regionstart, $regionend) = ($refhit =~ /([^_]+)\_([^_]+)\_([^_]+)/) ? ($1, $2, $3) : ('NA', 'NA', 'NA');

        push @{$mash_hits{$readname}}, { 'readlength' => $readlength,
                                         'refhit' => $refhit,
                                         'perciden' => $percident,
                                         'hitlength' => $readend - $readstart,
                                         'chrom' => $chrom,
                                         'regionstart' => $regionstart,
                                         'regionend' => $regionend } if ($readend - $readstart > 500);
    }
    close MASH;

    return {%mash_hits};
}
