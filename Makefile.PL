use ExtUtils::MakeMaker;

WriteMakefile(
    NAME               => 'check_for_gold_in_short_seqs.pl',
    VERSION_FROM       => 'check_for_gold_in_short_seqs.pl',
    LICENSE            => 'gpl_3',
    CONFIGURE_REQUIRES => { 'ExtUtils::MakeMaker' => '0', },
    PREREQ_PM          => {'JFR::Fasta => 1.0'}, 
    EXE_FILES          => ['check_for_gold_in_short_seqs.pl']);

