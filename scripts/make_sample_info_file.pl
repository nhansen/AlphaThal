#!/usr/local/bin/perl -w
#

use File::Path qw( make_path );
use Getopt::Long;

use strict;

our $SRAFILE = '/data/nhansen/HERV_K_catalog/SRAHumanPacBioWGSRuns.linux.051119.txt';
our $GENOMESIZE=3095693983;

my $Usage = qq!make_sample_info_file.pl <refgenotype directory>\n!;

$#ARGV==0
    or die "$Usage";

my $refgenotypedir = $ARGV[0];

my $rh_sample_info = retrieve_sample_accessions();

opendir SAMPLES, $refgenotypedir
    or die "Couldn\'t open $refgenotypedir for reading: $!\n";
my @sampledirs = readdir SAMPLES;
closedir SAMPLES;

foreach my $biosample (keys %{$rh_sample_info}) {
    my ($dirname) = grep /\.$biosample$/, @sampledirs;
    next if !($dirname);
    my $coverage = $rh_sample_info->{$biosample}->{total_bases}/$GENOMESIZE;
    $coverage = sprintf("%4.3f", $coverage);
    my $platforms = join ';', @{$rh_sample_info->{$biosample}->{platforms}};
    my $release_dates = join ';', sort {my $ayear=($a=~/\/(\d+)/) ? $1 : 0; my $byear=($b=~/\/(\d+)/) ? $1 : 0; return $ayear<=>$byear} @{$rh_sample_info->{$biosample}->{release_dates}};
    my $sample = ($dirname =~ /^(\S+)\.$biosample$/) ? $1 : $rh_sample_info->{$biosample}->{sample};
    print "$biosample\t$sample\t$release_dates\t$platforms\t$coverage\n";
}

sub retrieve_sample_accessions {

    open SRA, $SRAFILE
        or die "Couldn\'t open $SRAFILE: $!\n";

    my %sample_info = ();
    my @header_fields;
    while (<SRA>) {
        chomp;
        if (/^Run/) {
            @header_fields = split /\t/, $_;
        }
        else {
            my @srr_fields = split /\t/, $_;
            my %this_srr = ();
            for (my $i=0; $i<=$#header_fields; $i++) {
                $this_srr{$header_fields[$i]} = $srr_fields[$i];
            }
            my $release_date = $this_srr{'ReleaseDate'};
            $release_date =~ s/\s.*//;
            if ($release_date =~ m:(\d+)/(\d+)/(\d+):) {
                $release_date = "$1/$3";
            }
            my $biosample = $this_srr{'BioSample'};
            my $sample = $this_srr{'SampleName'};
            my $bases = $this_srr{'bases'};
            my $platform = $this_srr{'Platform'};
            my $model = $this_srr{'model'};
            $model = ($model) ? " $model" : "";
            $sample_info{$biosample}->{sample} = $sample;
            if (!grep {$_ eq $release_date} @{$sample_info{$biosample}->{release_dates}}) {
                push @{$sample_info{$biosample}->{release_dates}}, $release_date;
            }
            if (!grep {$_ eq "$platform $model"} @{$sample_info{$biosample}->{platforms}}) {
                push @{$sample_info{$biosample}->{platforms}}, "$platform $model";
            }
            $sample_info{$biosample}->{total_bases} += $bases;
        }
    }
    close SRA;

    return {%sample_info};
}
