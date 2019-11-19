#!/usr/local/bin/perl -w
#

use GTB::File qw(Open);
use Getopt::Long;

use strict;

our $PROJECTDIR='/data/nhansen/HERV_K_catalog/discovery2';
our $SCRIPTSDIR='/data/nhansen/HERV_K_catalog/discovery2/scripts';
our $MERGEDDIR='/data/nhansen/HERV_K_catalog/discovery2/refgenotype/mergedSVs';

my $Usage = qq!assemble_alternate_alleles.pl <hervk id> <merge directory>\n!;

my $minaltreads = 6;
my $maxmemory = '6g';
my $slurmthreads = 2;
GetOptions("minaltreads=i" => \$minaltreads, "maxmemory=s" => \$maxmemory, "slurmthreads=i" => \$slurmthreads );

$#ARGV==1
    or die "$Usage";

my $hervk_id = $ARGV[0];
my $mergedir = $ARGV[1];

my $top_dir = "$mergedir";

my $unionvcf = "$mergedir/merge/$hervk_id.merged.vcf";
my $rh_sv_readnames = read_union_vcf($unionvcf);
my $merged_file = "$mergedir/merge/$hervk_id.merged.sorted.clustered.svs.tdf";
my $rh_alt_alleles = read_hervk_altalleles($hervk_id, $top_dir, $merged_file, $rh_sv_readnames);

my $num_alleles = keys %{$rh_alt_alleles};

if (!$num_alleles) { # no alleles to process
    exit 0;
}

my $alleleplural = ($num_alleles == 1) ? '' : 's';
print "$hervk_id has $num_alleles allele$alleleplural!\n";

########END MAIN#######

sub read_union_vcf {
    my $union_vcf = shift;

    my %readname = ();
    my $vcf_fh = Open($union_vcf);
    while (<$vcf_fh>) {
        next if (/^#/);
        my ($id, $pos, $sv_id, $ref, $alt, $qual, $filter, $info, $remainder) = split /\t/, $_;
        my @info_strings = split /;/, $info;
        my ($altfield) = grep /^ALTPOS=/, @info_strings;
        if ($altfield =~ /ALTPOS=([^:]+):/) {
            my $readname = $1;
            $readname{$sv_id} = $readname;
        }
    }
    close $vcf_fh;
    return {%readname};

}

sub read_hervk_altalleles {
    my $hervk_id = shift;
    my $top_dir = shift;
    my $merged_file = shift;
    my $rh_sv_readnames = shift;

    my %allele_info = ();
    my $altid_number = 1;
    open ALTS, $merged_file
        or die "Couldn\'t open $merged_file for reading: $!\n";
    while (<ALTS>) {
        next if (/^#/);
        chomp;
        my ($hid, $pos, $svid, $reflength, $altlength,
            $qual, $end, $svtype, $reptype, $svlen,
            $bslength, $refwidened, $altwidened, $altpos,
            $clusterids, $numclustersvs, $exactids,
            $numexactsvs, $maxshiftdist, $maxsizediff,
            $maxeditdist) = split /\t/, $_;

        my @cluster_ids = split /:/, $clusterids;

        # split cluster corrected reads into their respective samples
        my $rh_clustersvs = {};
        foreach my $clustersv (@cluster_ids) {
            my ($sample, $svname) = ($clustersv =~ /^(\S+\.SAMN\d+)\.(\S+)$/) ? ($1, $2) : ('NA', 'NA');
            if ($sample eq 'NA') {
                print STDERR "Skipping unrecognized cluster member name $clustersv!\n";
                next;
            }
            push @{$rh_clustersvs->{$sample}}, $svname;
        }

        my %thisallele_data = ();
        foreach my $sample (keys %{$rh_clustersvs}) {
            my $num_reads = @{$rh_clustersvs->{$sample}};
            next if ($num_reads < $minaltreads);
            print "$sample has $num_reads alt reads for $hervk_id.ALT$altid_number\n";
            $thisallele_data{$sample} = [@{$rh_clustersvs->{$sample}}];
        }

        my $num_samples = grep !/HG03125/, keys %thisallele_data;
        if ($num_samples) {
            $allele_info{"$hervk_id.ALT$altid_number"} = {%thisallele_data};

            my $ctg_fasta = run_wtdbg($top_dir, $hervk_id, $altid_number, {%thisallele_data}, $rh_sv_readnames);

            if (-e $ctg_fasta) {
                system("samtools faidx $ctg_fasta");
                my $no_contigs = `grep -c '>' $ctg_fasta`;
                if ($no_contigs == 1) {
                    print STDERR "Successfully assembled one contig\n";
                }
                else {
                    print STDERR "FAILED: assembled $no_contigs contigs\n";
                }
            }
            else {
                print STDERR "FAILED: no assembled contig file\n";
            }

            $altid_number++;
        }
    }
    close ALTS;

    return {%allele_info};
}

sub run_wtdbg {
    my $top_dir = shift;
    my $hervk_id = shift;
    my $altid_number = shift;
    my $rh_allele_data = shift;
    my $rh_sv_readnames = shift;

    mkdir("$top_dir/assemble/$hervk_id.ALT$altid_number");
    my $hervfasta = "$top_dir/assemble/$hervk_id.ALT$altid_number/$hervk_id.ALT$altid_number.reads.fasta";

    my $hervfasta_fh = Open($hervfasta, "w");

    foreach my $sample (keys %{$rh_allele_data}) {
        my @cluster_svs = @{$rh_allele_data->{$sample}};
        my $correctedreadfasta = "$PROJECTDIR/refgenotype/$sample/canu_correct/$hervk_id/$hervk_id.correctedReads.fasta";
        foreach my $sv (@cluster_svs) {
            print "Retrieving readname for sv $sample.$sv\n";
            my $read = $rh_sv_readnames->{"$sample.$sv"};
            my $readfasta = `samtools faidx $correctedreadfasta $read`;
            foreach my $line (split /\n/, $readfasta) {
                if ($line =~ />\s*(\S+)/) {
                    my $readname = "$sample.$1";
                    print $hervfasta_fh ">$readname\n";
                }
                else {
                    print $hervfasta_fh "$line\n";
                }
            }
        }
    }
    close $hervfasta_fh;

    chdir("$top_dir/assemble/$hervk_id.ALT$altid_number");
    my $cmd = "wtdbg2 -i $hervfasta -fo $hervk_id.ALT$altid_number -g 40000 -x rs -t 2";
    print "$cmd\n";
    system($cmd);
    $cmd = "wtpoa-cns -t 2 -i $hervk_id.ALT$altid_number.ctg.lay.gz -fo $hervk_id.ALT$altid_number.ctg.fa";
    print "$cmd\n";
    system($cmd);

    return "$top_dir/assemble/$hervk_id.ALT$altid_number/$hervk_id.ALT$altid_number.ctg.fa";
}

sub run_canu {
    my $top_dir = shift;
    my $hervk_id = shift;
    my $altid_number = shift;
    my $rh_allele_data = shift;
    my $rh_sv_readnames = shift;

    mkdir("$top_dir/canu/$hervk_id.ALT$altid_number");
    my $hervfasta = "$top_dir/canu/$hervk_id.ALT$altid_number/$hervk_id.ALT$altid_number.reads.fasta";

    foreach my $sample (keys %{$rh_allele_data}) {
        my @cluster_svs = @{$rh_allele_data->{$sample}};
        my $correctedreadfasta = "$PROJECTDIR/refgenotype/$sample/canu_correct/$hervk_id/$hervk_id.correctedReads.fasta";
        foreach my $sv (@cluster_svs) {
            print "Retrieving readname for sv $sample.$sv\n";
            my $read = $rh_sv_readnames->{"$sample.$sv"};
            system("samtools faidx $correctedreadfasta $read >> $hervfasta");
        }
    }

    my $cmd = "canu -genomeSize=40k -assemble -pacbio-raw $hervfasta usegrid=0 corOutCoverage=100 -maxMemory=$maxmemory -maxThreads=$slurmthreads -p $hervk_id.ALT$altid_number -d $top_dir/canu/$hervk_id.ALT$altid_number";
    print "$cmd\n";
    system($cmd);

}

