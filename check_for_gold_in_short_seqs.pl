#!/usr/bin/perl

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

use strict;
use warnings;
use JFR::Fasta;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use Data::Dumper;

our $PROGRAM_NAME = 'check_for_gold_in_short_seqs.pl';
our $VERSION = 0.01;
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';
our $BLAT = 'pblat';
our $DEFAULT_STEPSIZE = 5;
our $DEFAULT_MINIDENTITY = 93;

# COMMENT BEFORE DISTRIBUTING
#use Storable qw(store retrieve);
#our $TMP_STORE_LONG = 'storable.rh_long_blat';
#store($rh_long_blat,$TMP_STORE_LONG);
#my $rh_long_blat = retrieve($TMP_STORE_LONG);
# END COMMENT

MAIN: {
    my $cmd_line = join " ", $0, @ARGV;
    my $rh_opts = process_options();        
    my $log_fh  = IO::File->new(">$rh_opts->{'out'}.log");
    print_cmd_and_version($cmd_line,$log_fh);
    my $blat_cmd = get_blat_cmd($BLAT,$rh_opts);
    my $long_blat_out = run_long_blat($blat_cmd,$rh_opts,$log_fh);
    my $rh_long_blat = get_blat($long_blat_out);
    my $masked = mask_transcripts($rh_opts,$rh_long_blat);
    my $short_blat = run_masked_transcriptome_blat($blat_cmd,$rh_opts,$masked,$log_fh);
    my $rh_results = get_results($short_blat);
    write_fasta($rh_opts,$rh_results,$log_fh);
    write_pairs_file($rh_opts->{'out'},$rh_results->{'pairs'},$log_fh);
    print_closing_msg($log_fh,$rh_opts);
}

sub print_cmd_and_version {
    my $commandline = shift;
    my $log_fh = shift;
    print $log_fh "$PROGRAM_NAME version $VERSION\n\n";
    print $log_fh "CMD: $commandline\n\n";
}

sub print_closing_msg {
    my $log_fh = shift;
    my $rh_o   = shift;
    print $log_fh qq~If the output files above are empty check the following BLAT error files: 
    $rh_o->{'out'}.shortmasked.blat.err, 
    $rh_o->{'out'}.long.blat.err\n\n$PROGRAM_NAME completed.\n~; 
 }

sub write_pairs_file {
    my $pre      = shift;
    my $ra_pairs = shift; 
    my $log_fh   = shift;
    open OUT, ">$pre.pairs.csv" or die "cannot open $pre.pairs.csv:$!";
    foreach my $ra_p (@{$ra_pairs}) {
        print OUT "$ra_p->[0],$ra_p->[1]\n";
    }
    print $log_fh "CSV file mapping transcripts to short seqs:\n    $pre.pairs.csv\n\n";
}

sub write_fasta {
    my $rh_o       = shift;
    my $rh_results = shift;
    my $log_fh     = shift;
    my $short_out = "$rh_o->{'out'}.short.w_uniq_transcriptome_matches.fa";
    open OUT, ">$short_out" or die "cannot open >$short_out:$!";
    my $fp = JFR::Fasta->new($rh_o->{'short_fasta'});
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        next unless ($rh_results->{'short'}->{$id});
        print OUT "$rec->{'def'}\n$rec->{'seq'}\n";
    }
    close OUT;
    print $log_fh "short sequences that are mapped to by uniq transcript stretches:\n    $short_out\n\n";
    my $tran_out = "$rh_o->{'out'}.transcripts_aln_to_short.fa";
    open OUT2, ">$tran_out" or die "cannot open >$tran_out:$!";
    my $fp2 = JFR::Fasta->new($rh_o->{'transcriptome'});
    while (my $rec = $fp2->get_record()) {
        my $id2 = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        next unless ($rh_results->{'transcripts'}->{$id2});
        print OUT2 "$rec->{'def'}\n$rec->{'seq'}\n";
    }
    close OUT2;
    print $log_fh "transcripts that have stretches that map uniquely to short sequences:\n    $tran_out\n\n";
}  

sub get_results {
    my $blat = shift;
    my %results = ();
    open IN, $blat or die "cannot open $blat:$!";
    while (my $line = <IN>) {
        my @fields = split /\t/, $line;
        $results{'transcripts'}->{$fields[9]}++;
        $results{'short'}->{$fields[13]}++;
        push @{$results{'pairs'}}, [$fields[9],$fields[13]];
    }
    return \%results;
}

sub run_masked_transcriptome_blat {
    my $blat_cmd = shift;
    my $rh_o = shift;
    my $masked = shift;
    my $log_fh = shift;
    $blat_cmd .= " $rh_o->{'short_fasta'} $masked $rh_o->{'out'}.shortmasked.psl > $rh_o->{'out'}.shortmasked.blat.out 2> $rh_o->{'out'}.shortmasked.blat.err";
    system $blat_cmd;
    print $log_fh "BLAT_CMD: $blat_cmd\n\n";
    return "$rh_o->{'out'}.shortmasked.psl";
}

sub get_blat {
    my $file = shift;
    my %data = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next if ($line =~ m/^#/);
        my @fields = split/\t/, $line;
        push @{$data{$fields[9]}->{$fields[13]}}, [$fields[11],$fields[12]];
    }
    return \%data;
}

sub mask_transcripts {
    my $rh_o = shift;
    my $rh_blat = shift;
    my $masked_out = "$rh_o->{'out'}.masked_transcriptome.fa";
    open OUT, ">$masked_out" or die "cannot open >$masked_out:$!";
    my $fp = JFR::Fasta->new($rh_o->{'transcriptome'});
    while (my $rec = $fp->get_record()) {
        print OUT "$rec->{'def'}\n";
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        if ($rh_blat->{$id}) {
            my @seq = split /|/, $rec->{'seq'};
            my $ra_blat  = create_bit_array($rh_blat->{$id});
            for (my $i = 0; $i < @seq; $i++) {
                if ($ra_blat->[$i]) {
                    print OUT "N";
                } else {
                    print OUT "$seq[$i]";
                }
            }
            print OUT "\n";
        } else {
            print OUT "$rec->{'seq'}\n";
        }
    }
    return $masked_out;
}

sub create_bit_array {
    my $rh_data = shift;
    my @array = ();
    foreach my $key (keys %{$rh_data}) {
        foreach my $ra_c (@{$rh_data->{$key}}) {
            die "unexpected" if ($ra_c->[0] > $ra_c->[1]);
            for (my $i=$ra_c->[0]; $i <= $ra_c->[1]; $i++) {
                $array[$i]++;
            }
        }
    }
    return \@array;
}


sub run_long_blat {
    my $blat_cmd = shift;
    my $rh_o = shift;
    my $log_fh = shift;
    $blat_cmd .= " $rh_o->{'long_fasta'} $rh_o->{'transcriptome'} $rh_o->{'out'}.long.psl > $rh_o->{'out'}.long.blat.out 2> $rh_o->{'out'}.long.blat.err";
    system $blat_cmd;
    print $log_fh "BLAT_CMD: $blat_cmd\n\n";
    return "$rh_o->{'out'}.long.psl";
}

sub get_blat_cmd {
    my $blat = shift;
    my $rh_o = shift;
    $blat = $rh_o->{'useblat'} if ($rh_o->{'useblat'});
    $blat .= " -threads=$rh_o->{'threads'}" if $rh_o->{'threads'};
    if ($rh_o->{'blatoptions'}) {
        $blat .= " $rh_o->{'blatoptions'}";
    } else {
        $blat .= " -noHead";
        $blat .= " -stepSize=$DEFAULT_STEPSIZE";
        $blat .= " -minIdentity=$DEFAULT_MINIDENTITY";
    }
    return $blat;
}

sub process_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                              "version" => \$rh_opts->{'version'},
                         "long_fasta=s" => \$rh_opts->{'long_fasta'},
                        "short_fasta=s" => \$rh_opts->{'short_fasta'},
                      "transcriptome=s" => \$rh_opts->{'transcriptome'},
                             "use_blat" => \$rh_opts->{'use_blat'},
                       "blat_options=s" => \$rh_opts->{'blat_options'},
                            "threads=i" => \$rh_opts->{'threads'},
                                "out=s" => \$rh_opts->{'out'},
                                 "help" => \$rh_opts->{'help'});
    die "$VERSION\n" if ($rh_opts->{'version'});
    pod2usage({-exitval => 0, -verbose => 2}) if $rh_opts->{'help'};
    unless ($rh_opts->{'long_fasta'} && $rh_opts->{'short_fasta'} &&
            $rh_opts->{'out'}) {
        warn "  missing --short_fasta\n" unless ($rh_opts->{'short_fasta'});
        warn "  missing --long_fasta\n" unless ($rh_opts->{'long_fasta'});
        warn "  missing --transcriptome\n" unless ($rh_opts->{'transcriptome'});
        warn "  missing --out\n" unless ($rh_opts->{'out'});
        warn "\n";
        usage();
    }
    # check files are readable and non-zero
    unless (-r $rh_opts->{'long_fasta'} && -s $rh_opts->{'long_fasta'}) {
      die "test if $rh_opts->{'long_fasta'} (arg to --long_fasta) is readable and non-zero failed";
    }
    unless (-r $rh_opts->{'short_fasta'} && -s $rh_opts->{'short_fasta'}) {
      die "test if $rh_opts->{'short_fasta'} (arg to --short_fasta) is readable and non-zero failed";
    }
    unless (-r $rh_opts->{'transcriptome'} && -s $rh_opts->{'transcriptome'}) {
      die "test if $rh_opts->{'transcriptome'} (arg to --transcriptome) is readable and non-zero failed";
    }
    return $rh_opts;
}

sub usage {
    print qq~usage: check_for_gold_in_short_seqs.pl
    --short_fasta=FASTA_W_SHORT_CONTIGS
    --long_fasta=FASTA_OF_ASSEMBLY_W_SHORT_CONTIGS_REMOVED
    --transcriptome=FASTA_OF_TRANSCRIPTOME_FROM_SAME_SPECIES
    --out=PREFIX_USED_FOR_OUTFILES
   [--use_blat]
   [--blat_options=OPTIONS_TO_PASS_TO_BLAT]
   [--threads=NUMBER_OF_THREADS_FOR_PBLAT_IF_USING_PBLAT]
   [--help]
   [--version]\n~;
    exit;
}

__END__

=head1 NAME

B<check_for_gold_in_short_seqs.pl> - ID short contigs w/transcript alns

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

check_for_gold_in_short_seqs.pl --short_fasta=FASTA_W_SHORT_CONTIGS --long_fasta=FASTA_OF_ASSEMBLY_W_SHORT_CONTIGS_REMOVED --transcriptome=FASTA_OF_TRANSCRIPTOME_FROM_SAME_SPECIES --out=PREFIX_USED_FOR_OUTFILES [--use_blat] [--blat_options=OPTIONS_TO_PASS_TO_BLAT] [--threads=NUMBER_OF_THREADS_FOR_PBLAT_IF_USING_PBLAT] [--help] [--version]

=head1 DESCRIPTION

This program uses blat or (preferably) pblat to align transcriptome sequences to a genome assembly. It then masks the transcriptome and reblats it vs. any short sequences that were removed from the assembly (for whatever reason). The program produces a log file (the prefix you provide with --out followed by .log), which explains output.

=head1 BUGS

Please report them to <joseph.ryan@whitney.ufl.edu>

=head1 OPTIONS

=over 2

=item B<--short_fasta>

The FASTA file with short sequences from an assembly

=item B<--long_fasta>

The FASTA file with the long contigs/scaffolds of the main assembly

=item B<--out>

A prefix to use for the output FASTQ (or FASTA) files

=back 

=head1 COPYRIGHT

Copyright (C) 2020 Joseph F. Ryan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
