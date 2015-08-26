# PGRN_JuncBASE_postprocessing
JuncBASE post-processing code used in PGRN RNA-seq manuscript Chhibber, et al. 2015
(for use with JuncBASE v0.6)

-Run JuncBASE through the 'create alt splicing tables of read counts' step (createAS_CountTables.py or run_createAS_CountTables.py/combine_createAS_CountTables_by_chr.py).


-Group together related JuncBASE events ('juncEvents'). Uses code from JuncBASE v0.6 makeNonRedundantAS.py. Works in output of previous step. For time, 'intron_retention', 'jcn_only_AD', and 'jcn_only_AA' events were filtered out before running this step.

```
python JuncBASE_postprocessing.py 
Usage: JuncBASE_postprocessing.py [options]
Options:
  -h, --help            show this help message and exit
  --directory=DIRECTORY
                        directory of files to fix
  --prefix=PREFIX       prefix of files to fix
```

Example command from PGRN analysis:
```
python JuncBASE_postprocessing.py --directory juncBASE_subsampled/20M/ --prefix allsamples_filNoIR_
```
  
-Define mutually exclusive 'junction sets' for each group from 'juncEvents', and report length normalized read counts. For example, a simple cassette exon events is reported on one line in the default JuncBASE output (juncEvent), with read counts for each sample reported as X;x for exclusion;inclusion. This script transforms that into two lines, one for the inclusion of the exon and one for the exclusion of the exon, each with only their own read counts. Also, makes a bed file for easy viewing of groups and junction sets. Filter out complex events where it is not possible to get accurate read counts.

```
python JuncBASE_postprocessing_2.py 
Usage: JuncBASE_postprocessing_2.py [options]
Options:
  -h, --help            show this help message and exit
  -i INPUT              input file
  -e ENTROPY            entropy_scores.txt with K/N
  --e_alt               -e input is list of known introns (     )
  -o OUTPUT             output file
  -b BED                output bed file
  -c COMPLEX            output file of complex events
  -n NAME               track name (opt)
  --jcn_seq_len=JCN_SEQ_LEN
                        jcn_seq_len
```

Example command from PGRN analysis:

```
python JuncBASE_postprocessing_2.py -i juncBASE_subsampled/20M/allsamples_filNoIR_grouped_AS_exclusion_inclusion_counts_lenNorm.txt -e gencode.v12.AND_ensGene_130510_exonsOnly.combined_intronCoords.txt --e_alt -o juncBASE_subsampled/20M/allsamples_filNoIR_grouped_junctionSets_counts_lenNorm.txt -b juncBASE_subsampled/20M/allsamples_filNoIR_grouped_junctionSets.bed -c juncBASE_subsampled/20M/allsamples_filNoIR_grouped_AS_exclusion_inclusion_counts_lenNorm_complexEvents.txt -n allsamples_filNoIR_simple_subsampled20M
```

-Calculate 'percent spliced in' (PSI) values by running JuncBASE v0.6 compareSampleSets.py with the 'as_only' option (no differential testing done) on the output of Step 2 above.


-Determine the PSI values for each junctions set using those calculated in Step 4 above ('a' option) and the junctions sets defined in Step 3 above ('j' option).

```
python JuncBASE_postprocessing_3.py
Usage: JuncBASE_postprocessing_3.py [options]
Options:
  -h, --help       show this help message and exit
  -a ALL_PSI_OUT   all_psi_out
  -j JUNCTIONSETS  junctionSets
```

Example command from PGRN analysis:
```
python JuncBASE_postprocessing_3.py -a juncBASE_subsampled/20M/allsamples_filNoIR_grouped_all_psi_output.txt -j juncBASE_subsampled/20M/allsamples_filNoIR_grouped_junctionSets_counts_lenNorm.txt >juncBASE_subsampled/20M/allsamples_filNoIR_grouped_junctionSets_PSIs.txt
```


-Run JuncBASE v0.6 compareSampleSets.py on the output of Step 2 above for each pair of sample sets as desired. 


-Determine the p-value for differential splicing for each junctions set using those calculated in Step 6 above ('a' option) and the junctions sets defined in Step 3 above ('j' option).

```
python JuncBASE_postprocessing_4.py
Usage: JuncBASE_postprocessing_4.py [options]
Options:
  -h, --help       show this help message and exit
  -a ALL_PSI_OUT   dir of inputs
  -j JUNCTIONSETS  junctionSets
```

Example command from PGRN analysis:

```
python JuncBASE_postprocessing_4.py -a juncBASE_subsampled/20M/betweenTissue/ -j juncBASE_subsampled/20M/allsamples_filNoIR_grouped_junctionSets_counts_lenNorm.txt >juncBASE_subsampled/20M/allsamples_filNoIR_grouped_junctionSets_betweenTissue_pvals.txt
```
