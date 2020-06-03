### To run the sample files from this directory:
perl ../check_for_gold_in_short_seqs.pl \
    --short_fasta=short.fa
    --long_fasta=long.fa
    --transcriptome=transcriptome.fa
    --out=sample

## File and sequence explanations below

#### short.fa
#### seqs that transcriptomes align to or best BLAST hit if no alignment
Bova5.0.lt200.001669
Bova5.0.lt200.003395
Bova5.0.lt200.116202
Bova5.0.lt200.156701
Bova5.0.lt200.008933 # has a transcribed stretch not in long.fa

#### long.fa
#### seqs that transcriptomes align to or best BLAST hit if no alignment
Bova5.00789
Bova5.01220
Bova5.02079
Bova5.03464
Bova5.35507
Bova5.12533

#### transcriptome
Bero_ovat_ALL.000001 = spread across two long scaffolds
Bero_ovat_ALL.000002 = encompassed in one long scaffold
Bero_ovat_ALL.000003 = not present in long or short scaffolds
Bero_ovat_ALL.000017 = partially in long, partially in short



