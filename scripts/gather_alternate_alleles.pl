#!/usr/local/bin/perl -w
#
# Script to gather merged SV calls from corrected reads, and construct phased
# alternate allele haplotypes. Finally, align uncorrected reads to both 
# alternate and reference haplotypes.

use Getopt::Long;

use strict;

our $LONGREADTOPDIR=$ENV{'LONGREADTOPDIR'} || '/data/nhansen/HERV_K_catalog/SVA_discovery';
our $MERGEDDIR=$LONGREADTOPDIR.'/refgenotype/mergedSVs';
our $REFGENODIR=$LONGREADTOPDIR.'/refgenotype';

my $Usage = qq!gather_alternate_alleles.pl <element id> <merge directory name>\n!;

my $minaltreads = 5;

GetOptions("minaltreads=i" => \$minaltreads );

$#ARGV==1
    or die "$Usage";

my $element_id = $ARGV[0];
my $mergedirname = $ARGV[1];

# read in structural variants with minimum number of reads:

my ($ra_all_svs, $rh_hap_reads) = read_altalleles($element_id, "$MERGEDDIR/$mergedirname/targetvcfs/$element_id.merged.sorted.clustered.vcf.gz");

my $hap_info_file = "$MERGEDDIR/$mergedirname/altalleleinfo/$element_id.althapreads.txt";
my $haplotype_fasta_file = "$MERGEDDIR/$mergedirname/altalleleinfo/$element_id.althaps.fasta";
my $ref_fasta = "$LONGREADTOPDIR/prep/allele_sequences/$element_id.widerregion.fasta";
my $correctedreads_fasta_file = "$MERGEDDIR/$mergedirname/altalleleinfo/$element_id.correctedreads.fasta";

my $ref_allele = read_ref_allele($ref_fasta);

construct_alternate_haplotypes($ra_all_svs, $rh_hap_reads, $hap_info_file, $haplotype_fasta_file, $ref_allele);

write_corrected_read_fasta($rh_hap_reads, $correctedreads_fasta_file, $REFGENODIR);

########END MAIN#######

sub read_altalleles {
    my $element_id = shift;
    my $merged_vcf_file = shift;

    my $rh_sv_reads = {};
    my %all_reads = ();
    my @all_svs = (); # ordered list of all SVs for this targeted element

    my $merged_vcf_string = "gunzip -c $merged_vcf_file | ";
    open ALTS, $merged_vcf_string
        or die "Couldn\'t open $merged_vcf_string for reading: $!\n";
    while (<ALTS>) {
        next if (/^#/);
        chomp;

        my ($chrom, $pos, $svid, $refstring, $altstring,
            $qual, $filter, $info, $rest) = split /\t/, $_;

        my @info_fields = split /;/, $info;
        my @numberfield = grep /NumClusterSVs=\d+/, @info_fields;
        my @clusterfields = grep /ClusterIDs=.*/, @info_fields;

        my $numclustersvs = ($numberfield[0] =~ /(\d+)/) ? $1 : 0;
        $clusterfields[0] =~ s/ClusterIDs=//;
        my @cluster_reads = split /:/, $clusterfields[0];
        my $numclusterreads = @cluster_reads;
        if ($numclustersvs != $numclusterreads) {
            die "Improper VCF file: $numclustersvs does not equal $numclusterreads for SV $svid!\n";
        }

        next if ($minaltreads && $numclustersvs < $minaltreads);

        map { s/\.\d+$// } @cluster_reads;

        $rh_sv_reads->{$svid} = [@cluster_reads];
        push @all_svs, {'svid' => $svid, 'pos' => $pos, 'ref' => $refstring, 'alt' => $altstring};
        foreach my $read (@cluster_reads) {
            $all_reads{$read} = 1;
        }
    }

    close ALTS;

    # phase reads, finding different haplotypes:
    # sort by position first:
 
    @all_svs = sort {$a->{pos} <=> $b->{pos}} @all_svs;

    my %haplotype_reads = ();
    foreach my $read (keys %all_reads) {
        my $hapstring = ''; # will append 0 if ref, 1 if alt for each SV
        foreach my $rh_sv (@all_svs) {
            my $sv = $rh_sv->{'svid'};
            if (grep {$_ eq $read} @{$rh_sv_reads->{$sv}}) {
                $hapstring .= '1';
            }
            else {
                $hapstring .= '0';
            }
        }
        # record hapstring:
        push @{$haplotype_reads{$hapstring}}, $read;
    }

    # write out results:
    
    foreach my $hapstring (keys %haplotype_reads) {
        print "$hapstring:\n";
        foreach my $read (@{$haplotype_reads{$hapstring}}) {
            print "\t$read\n";
        }
    }

    return ([@all_svs], {%haplotype_reads});
}

sub construct_alternate_haplotypes {
    my $ra_all_svs = shift;
    my $rh_hap_reads = shift;
    my $info_file = shift;
    my $fasta_file = shift;
    my $ref_allele = shift;

    open INFO, ">$info_file"
        or die "Couldn\'t open $info_file for writing: $!\n";
    open FASTA, ">$fasta_file"
        or die "Couldn\'t open $fasta_file for writing: $!\n";
    # Write reference allele:
    my $splitallele = $ref_allele;
    $splitallele =~ s/(.{50})/$1\n/g;
    chomp $splitallele;
    print FASTA ">$element_id.REF\n$splitallele\n";
    my $sv_number = 1;
    foreach my $hapstring (keys %{$rh_hap_reads}) {
        my $alt_allele = "$element_id.ALT$sv_number";
        my $reversehapstring = reverse $hapstring;
        my $allele = $ref_allele;
        my $index = 0;
        my @includedsvs = ();
        while (($reversehapstring ne '') && (my $nextgeno = chop $reversehapstring)) {
            if ($nextgeno eq '1') {
                $allele = apply_sv($allele, $ra_all_svs, $index);
                push @includedsvs, $ra_all_svs->[$index]->{'svid'};
                print INFO "HAP\t$alt_allele\t$ra_all_svs->[$index]->{'svid'}\n";
            }
            $index++;
        }
        $splitallele = $allele;
        $splitallele =~ s/(.{50})/$1\n/g;
        chomp $splitallele;
        print FASTA ">$alt_allele\n$splitallele\n";
        $sv_number++;
        my $svstring = join ':', @includedsvs;
        foreach my $read (@{$rh_hap_reads->{$hapstring}}) {
            print INFO "READ\t$alt_allele\t$read\n";
        }
    }
    close FASTA;
    close INFO;

}

sub write_corrected_read_fasta {
    my $rh_hap_reads = shift;
    my $fasta_file = shift;
    my $refgenodir = shift;

    opendir SAMPLES, $refgenodir
        or die "Couldn\'t open $refgenodir: $!\n";
    my @sampledirnames = grep { -e "$refgenodir/$_/canu_correct/$element_id/$element_id.correctedReads.fasta.gz" } readdir SAMPLES;

    closedir SAMPLES;

    open FASTA, ">$fasta_file"
        or die "Couldn\'t open $fasta_file for writing: $!\n";

    foreach my $sampledir (@sampledirnames) {
        open READS, "gunzip -c $refgenodir/$sampledir/canu_correct/$element_id/$element_id.correctedReads.fasta.gz | "
            or die "Couldn\'t open and gunzip $refgenodir/$sampledir/$element_id/$element_id.correctedReads.fasta.gz for reading: $!\n";

        while (<READS>) {
            if (/^>(\S+)/) {
                my $readname = $1;
                print FASTA ">$sampledir:$readname\n";
                my $readseq = <READS>;
                print FASTA "$readseq";
            }
        }
        close READS;
    }

    close FASTA;

}

sub apply_sv {
    my $allele = shift;
    my $ra_svs = shift;
    my $index = shift;

    my $pos = $ra_svs->[$index]->{pos};
    my $ref = $ra_svs->[$index]->{ref};
    my $reflength = length($ref);
    my $alt = $ra_svs->[$index]->{alt};
    my $altlength = length($alt);

    substr($allele, $pos-1, $reflength) = $alt;

    my $position_change = $altlength - $reflength;
    for (my $i=$index+1; $i<=$#{$ra_svs}; $i++) {
        $ra_svs->[$i]->{pos} += $position_change;
    }

    return $allele;
}

sub read_ref_allele {
    my $fasta = shift;

    open REF, $fasta
        or die "Couldn\'t open $fasta for reading: $!\n";

    my $ref_allele = '';
    while (<REF>) {
        next if (/>/);
        chomp;
        $ref_allele .= $_;
    }
    close REF;

    return $ref_allele;
}
