#!/usr/local/bin/perl -w
#
# Script to gather all necessary files for the long read viewer into 
# three directories: a "projectdata" directory, a "sampledata" 
# directory, and a "haplotypedata" directory.

use File::Path qw( make_path );
use Getopt::Long;

use strict;

our $LONGREADTOPDIR = $ENV{'LONGREADTOPDIR'};

my $Usage = qq!prepare_viewer_files.pl <merge directory name>\n!;

my $verbose = 0; # print out lots of info?

GetOptions( verbose => \$verbose );

$#ARGV==0
    or die "$Usage";

my $mergedirname = $ARGV[0];

if (!$LONGREADTOPDIR) {
    die "You must set the LONGREADTOPDIR environment variable before running this script!\n";
}

our $CURRENTDIR = $ENV{'PWD'};

if ($LONGREADTOPDIR ne $CURRENTDIR) {
    die "LONGREADTOPDIR variable is not set to current directory";
    print "LONGREADTOPDIR variable is not set to current directory--reset to current directory? ";
    my $answer = <STDIN>;
    chomp $answer;
    if ($answer =~ /^[Yy]/) {
        print "Resetting LONGREADTOPDIR to $CURRENTDIR!\n";
        $ENV{'LONGREADTOPDIR'} = $CURRENTDIR;
    }
    else {
        print "LONGREADTOPDIR will remain as $LONGREADTOPDIR!\n";
    }
}
our $MERGEDDIR=$LONGREADTOPDIR.'/refgenotype/mergedSVs';
our $REFGENODIR=$LONGREADTOPDIR.'/refgenotype';

# create necessary directories:
my $rh_dirs = make_directories($mergedirname);

# projectdata files:
create_projectdata_files($rh_dirs);

# sampledata bam links:
opendir SAMPLES, "$LONGREADTOPDIR/refgenotype"
    or die "Couldn\'t open directory $LONGREADTOPDIR/refgenotype for reading: $!\n";

my @sampledirs = grep { -e "$LONGREADTOPDIR/refgenotype/$_/allele_aligns/$_.genome.converted.sort.bam" } readdir SAMPLES;
closedir SAMPLES;

foreach my $sampledir (@sampledirs) {
    if (!(-e "$rh_dirs->{sampledatadir}/$sampledir.genome.converted.sort.bam")) {
        symlink("$LONGREADTOPDIR/refgenotype/$sampledir/allele_aligns/$sampledir.genome.converted.sort.bam", "$rh_dirs->{sampledatadir}/$sampledir.genome.converted.sort.bam");
    }
    if (!(-e "$rh_dirs->{sampledatadir}/$sampledir.genome.converted.sort.bam.bai")) {
        symlink("$LONGREADTOPDIR/refgenotype/$sampledir/allele_aligns/$sampledir.genome.converted.sort.bam.bai", "$rh_dirs->{sampledatadir}/$sampledir.genome.converted.sort.bam.bai");
    }
}

# haplotypedata files:
opendir HAPS, "$MERGEDDIR/$mergedirname/altalleleinfo"
    or die "Couldn\'t open directory $MERGEDDIR/$mergedirname/altalleleinfo for reading: $!\n";

my @althap_fastas = grep /\.althaps.fasta$/, readdir HAPS;
closedir HAPS;

open HAPINFO, ">$rh_dirs->{haplotypedatadir}/target_haplotype_info.txt"
    or die "Couldn\'t open file $rh_dirs->{haplotypedatadir}/target_haplotype_info.txt for writing: $!\n";

open SVVCF, ">$rh_dirs->{haplotypedatadir}/genomic_svs.vcf"
    or die "Couldn\'t open file $rh_dirs->{haplotypedatadir}/genomic_svs.vcf for writing: $!\n";

open SVINFO, ">$rh_dirs->{haplotypedatadir}/genomic_sv_info.txt"
    or die "Couldn\'t open file $rh_dirs->{haplotypedatadir}/genomic_sv_info.txt for writing: $!\n";

open ALLHAPS, "| gzip -c > $rh_dirs->{haplotypedatadir}/target_haplotypes.fasta.gz"
    or die "Couldn\'t open $rh_dirs->{haplotypedatadir}/target_haplotypes.fasta.gz for writing: $!\n";

open HAPCOUNTS, ">$rh_dirs->{haplotypedatadir}/haplotype_read_counts.txt"
    or die "Couldn\'t open $rh_dirs->{haplotypedatadir}/haplotype_read_counts.txt for writing: $!\n";

my $target_bed_file = "$LONGREADTOPDIR/prep/ref_and_nonref.wideregions.bed";
my %target_data = read_target_positions($target_bed_file);

foreach my $althapfasta (sort @althap_fastas) {
    open TARGET, "$MERGEDDIR/$mergedirname/altalleleinfo/$althapfasta"
        or die "Couldn\'t open $MERGEDDIR/$mergedirname/altalleleinfo/$althapfasta for reading: $!\n";
    while (<TARGET>) {
        print ALLHAPS $_;
    }
    close TARGET;

    my $hapinfo_file = $althapfasta;
    $hapinfo_file =~ s/althaps\.fasta$/althapreads.txt/;

    open HAPS, "$MERGEDDIR/$mergedirname/altalleleinfo/$hapinfo_file"
        or die "Couldn\'t open $MERGEDDIR/$mergedirname/altalleleinfo/$hapinfo_file for reading: $!\n";

    my %hap_svs = ();
    while (<HAPS>) {
        if (/^HAP/) {
            chomp;
            my @fields = split /\t/, $_;
            my $targetid = $fields[1];
            $targetid =~ s/\..*$//;
            print HAPINFO "$targetid\t$fields[1]\t$fields[2]\t$fields[3]\n";
            $hap_svs{$fields[2]} = $fields[3];
        }
    }
    close HAPS;

    # file to pull SVs from:
    my $element_id = ($althapfasta =~ /^([^.]+)\./) ? $1 : 'NA';
    if ($element_id eq 'NA') {
        die "Unable to determine element id from haplotype fasta name $althapfasta!\n";
    }
    my $merged_sv_vcf = "$MERGEDDIR/$mergedirname/targetvcfs/$element_id.merged.sorted.clustered.vcf";

    my $vcf_open_string = (-e $merged_sv_vcf) ? "$merged_sv_vcf" : "gunzip -c $merged_sv_vcf.gz | ";

    open VCF, $vcf_open_string
        or die "Couldn\'t open $vcf_open_string: $!\n";
    while (<VCF>) {
        next if (/^#/);
        chomp;
        my @fields = split /\t/, $_;
        if ($hap_svs{$fields[2]}) {
            $fields[2] = $hap_svs{$fields[2]};
            my $rh_element_data = $target_data{$fields[0]};
            $fields[0] = $rh_element_data->{chrom};
            $fields[1] += $rh_element_data->{start} - 1;
            my $new_end=0;
            if ($fields[7] =~ /END=(\d+);/) {
                my $old_end = $1;
                $new_end = $old_end + $rh_element_data->{start} - 1;
                $fields[7] =~ s/END=$old_end;/END=$new_end;/;
            }
            my $vcf_record = join "\t", @fields;
            print SVVCF "$vcf_record\n";
            my $type = ($fields[7] =~ /REPTYPE=([^;]+);/) ? $1 : 'NA';
            my $length = ($fields[7] =~ /SVLEN=\-{0,1}(\d+);/) ? $1 : 'NA';
            print SVINFO "$element_id\t$fields[2]\t$fields[0]\t$fields[1]\t$new_end\t$type\t$length\n";
        }
    }
    close VCF;

    # count reads for each sample aligning to each haplotype:
    my $cm_file = $althapfasta;
    $cm_file =~ s/\.althaps.fasta/.reads_vs_althaps.cm.out/;

    my %max_score = ();
    my %max_match = ();
    if (-e "$MERGEDDIR/$mergedirname/altalleleinfo/$cm_file") {
        open CM, "$MERGEDDIR/$mergedirname/altalleleinfo/$cm_file"
            or die "Couldn\'t open $MERGEDDIR/$mergedirname/altalleleinfo/$cm_file!";
        while (<CM>) {
            if (/^ALIGNMENT\s+(\d+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+C{0,1}\s*(\S+)/) {
                my ($score, $read, $hap) = ($1, $2, $3);
                my ($sample, $readname) = ($read =~ /(\S+):(\S+)/) ? ($1, $2) : ('NA', 'NA');
                if (!$max_score{$sample} || !$max_score{$sample}->{$readname} || 
                      $max_score{$sample}->{$readname} < $score) {
                    $max_score{$sample}->{$readname} = $score;
                    $max_match{$sample}->{$readname} = $hap;
                }
            }
        }
        close CM;
        foreach my $sample (keys %max_score) {
            my %hap_counts = ();
            foreach my $readname (keys %{$max_match{$sample}}) {
                my $hap = $max_match{$sample}->{$readname};
                print "$sample $readname matches $hap\n";
                $hap_counts{$hap}++;
            }
            my $hapstring = join '/', keys %hap_counts;
            my $countstring = join '/', values %hap_counts;
            print HAPCOUNTS "$element_id\t$sample\t$hapstring\t$countstring\n";
        }
    }
}
close ALLHAPS;
close HAPINFO;
close SVVCF;
close SVINFO;
close HAPCOUNTS;

########END MAIN#######

sub make_directories {
    my $mergedirname = shift;

    my %directories = ();

    my $mergedir = "$MERGEDDIR/$mergedirname";
    ($directories{'mergedir'}) = (-e "$mergedir") ? $mergedir :  make_path $mergedir;
   
    ($directories{'viewerdir'}) = (-e "$mergedir/viewerfiles") ? "$mergedir/viewerfiles" : make_path "$mergedir/viewerfiles";

    my $projectdatadir = "$mergedir/viewerfiles/projectdata"; 
    ($directories{'projectdatadir'}) = (-e "$projectdatadir") ? "$projectdatadir" : make_path "$projectdatadir";

    my $sampledir = "$mergedir/viewerfiles/sampledata"; 
    ($directories{'sampledatadir'}) = (-e "$sampledir") ? "$sampledir" : make_path "$sampledir";

    my $haplotypedir = "$mergedir/viewerfiles/haplotypedata"; 
    ($directories{'haplotypedatadir'}) = (-e "$haplotypedir") ? "$haplotypedir" : make_path "$haplotypedir";

    return {%directories};
}

sub create_projectdata_files {
    my $rh_dirs = shift;

    system("$LONGREADTOPDIR/scripts/make_sample_info_file.pl $LONGREADTOPDIR/refgenotype > $rh_dirs->{projectdatadir}/sample_info.txt");
    # tab-delimited file with: chrom, start, end, annotationtype (for coloring and display along ref in read view)
    symlink("$LONGREADTOPDIR/prep/target_annotations.bed", "$rh_dirs->{projectdatadir}/target_annotations.bed");

    # colors for annotation elements:
    my @annot_colors = ("red", "pink", "orange", "blue", "purple", "brown", "yellow");
    my $command = "awk '{print \$1}' $LONGREADTOPDIR/prep/target_annotations.bed | sort | uniq | ";
    open ANNOTS, $command
        or die "Couldn\'t open $command for execution: $!\n";

    my $annotcolorfile = "$rh_dirs->{projectdatadir}/annotation_colors.txt";
    open ANNOCOLOR, ">$annotcolorfile"
        or die "Couldn\'t open $annotcolorfile for writing: $!\n";

    my @remaining_colors = @annot_colors;
    while (<ANNOTS>) {
        chomp;
        my $annot = $_;
        if (!@remaining_colors) {
            @remaining_colors = @annot_colors;
        }
        my $color = shift @remaining_colors;
        print ANNOCOLOR "$annot\t$color\n";
    }
    close ANNOTS;
    close ANNOCOLOR;

    # tab-delimited file with: target name, alias
    symlink("$LONGREADTOPDIR/prep/target_aliases.txt", "$rh_dirs->{projectdatadir}/target_aliases.txt");
    symlink("$LONGREADTOPDIR/prep/ref_and_nonref_target_regions.withgenes.bed", "$rh_dirs->{projectdatadir}/ref_and_nonref_target_regions.withgenes.bed");

}

sub read_target_positions {
    my $bedfile = shift;

    my %target_data = ();
    open TARGETS, $bedfile
        or die "Couldn\'t open $bedfile for reading: $!\n";

    while (<TARGETS>) {
        chomp;
        my ($chrom, $start, $end, $elementid) = split /\t/, $_;
        $target_data{$elementid} = {'chrom' => $chrom, 'start' => $start, 'end' => $end};
    }
    close TARGETS;

    return %target_data;
}
