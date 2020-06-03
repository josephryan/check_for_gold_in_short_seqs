# check_for_gold_in_short_seqs
It is common practice before finalizing genome assemblies to remove short seqs. This script uses transcripts to identify bits of genes missing in the full assembly but present in the short. Note:GenBank requires removal of seqs shorter than 200.

## REQUIREMENTS

1. perl
2. JFR-PerlModules https://github.com/josephryan/JFR-PerlModules
3. blat or (preferably) pblat 
   (#if using blat instead of pblat need to use --use_blat option)

## INSTALL

To install type the following:

    perl Makefile.PL

    make

    make install

To install without root privelages try:

    perl Makefile.PL PREFIX=/home/myuser/scripts

    make

    make install

## USAGE

check_for_gold_in_short_seqs.pl \
    --short_fasta=FASTA_W_SHORT_CONTIGS \
    --long_fasta=FASTA_OF_ASSEMBLY_W_SHORT_CONTIGS_REMOVED \
    --transcriptome=FASTA_OF_TRANSCRIPTOME_FROM_SAME_SPECIES \
    --out=PREFIX_USED_FOR_OUTFILES \
   [--use_blat]
   [--blat_options=OPTIONS_TO_PASS_TO_BLAT]
   [--threads=NUMBER_OF_THREADS_FOR_PBLAT_IF_USING_PBLAT]
   [--help]
   [--version]   

## DOCUMENTATION

Run the following for documentation:
perldoc check_for_gold_in_short_seqs.pl

COPYRIGHT AND LICENCE
------------

Copyright (C) 2019 by Joseph Ryan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program in the file LICENSE.  If not, see
http://www.gnu.org/licenses/.

