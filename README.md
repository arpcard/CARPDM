# CARPDM

## The Comprehensive Antibiotic Resistance Probe Design Machine

This program uses BLASTN to design probesets for targeted enrichment of the input sequences from DNA sequencing libraries. If desired, it can also modify these sequences to be an oligo pool, which is used for the in-house synthesis of the probeset, vastly reducing cost. Finally, it outputs several analysis graphics and summary statistics to guage the effectivenes of the probe design strategy.

## License 

Use or reproduction of these materials, in whole or in part, by any commercial organization whether or not for non-commercial (including research) or commercial purposes is prohibited, except with written permission of McMaster University. Commercial uses are offered only pursuant to a written license and user fee. To obtain permission and begin the licensing process, see the [CARD website](https://card.mcmaster.ca/about)

## Citation

Hackenberger et al. 2024. CARPDM: cost-effective antibiotic resistome profiling of metagenomic samples using targeted enrichment. bioRxiv 2024.03.27.587061

## Support & Bug Reports

Please create a [github issue](https://github.com/arpcard/CARPDM/issues) for any problems with the program.

You can email the CARD curators or developers directly at card@mcmaster.ca.

*******************************************************************************

## Installation

### Dependencies via Conda

If not already installed, follow the documentation to [install Miniconda](https://docs.anaconda.com/miniconda/)

Clone the github repository or download the carpdm_env.yml file.

```console
git clone https://github.com/arpcard/CARPDM.git
```

Create a conda environment from the carpdm_env.yml file.

```console
conda env create -f CARPDM/carpdm_env.yml
```

Activate the conda environment before running the script

```console
conda activate carpdm
```

******************************************************************************

## Running the CARPDM script

If desired, add the directory containing the carpdm.py script to your path by modifying your .bashrc file.

Add execute permissions to the carpdm.py script

```console
chmod a+x carpdm.py
```

See the full help menu by running 

```console
carpdm.py --help

usage: carpdm.py [-h] -i INPUT_FASTA [-p PROBE_NUM_CUTOFF] [-s TILING_STEP] [-l PROBE_LENGTH]
                 [-m MELT_TEMP] [-b BASENAME] [-o OUTPUT_DIR] [-f FILTER_FASTA | -d FILTER_DB]
                 [-t NUM_THREADS] [-n] [-c] [-v]

Design probes against an input fasta according to design parameters, then pass through several
filters to remove off-target enrichment (nt filter) and probe redundancy (self filter)

options:
  -h, --help            show this help message and exit
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Path to the fasta that contains sequences against which probes are to be
                        designed. Required
  -p PROBE_NUM_CUTOFF, --probe_num_cutoff PROBE_NUM_CUTOFF
                        Maximum number of final probes allowed in the probeset. The redundancy filter
                        will iterate until it reaches below this cutoff. Default = 42000
  -s TILING_STEP, --tiling_step TILING_STEP
                        Number of bases by which to offset initial probe tiling. Higher values make
                        computation less intense, though may result in sparser coverage. Default = 4
  -l PROBE_LENGTH, --probe_length PROBE_LENGTH
                        Probe length to create. Note that if >80 (default), probes will be too long
                        for the base price of a Twist oligo pool after addition of amplification and
                        transcription primers for in-house synthesis. Default = 80
  -m MELT_TEMP, --melt_temp MELT_TEMP
                        Minimum melting temperature for probes to be considered during basic filter.
                        Default = 50
  -b BASENAME, --basename BASENAME
                        File basename
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory
  -f FILTER_FASTA, --filter_fasta FILTER_FASTA
                        Fasta file to filter against. This creates a blast database from the provided
                        file to filter probes against after the basic filter. Probes with over
                        >(probe_length * 0.625) identities against anything in this fasta file will
                        be removed
  -d FILTER_DB, --filter_db FILTER_DB
                        Premade blast database to filter against. Probes with >(probe_length * 0.625)
                        identities against anything in this fasta file will be removed, unless the
                        specified database is the nt database (basename == "nt"). In that case it
                        will only remove probes that align with that condition to non-bacterial
                        accessions, that also do not have >(probe_length * 0.975) identities
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads to be used during BLAST searching
  -n, --no_o_pool       
                        If included, omits steps to design template oligo pool. Only include if
                        ordering RNA probes directly from the manufacturer.
  -c, --clean           
                        If included, removes all files except for final baits, oligo pools, and analysis plots
  -v, --version         
                        show program's version number and exit
```

We recommend a test run using the clinically_relevant_amr.fna file. Optionally, you might include a fasta file or BLAST database against which to filter. When filtering against the *E. coli* K12 reference genome, this command took just over 30 minutes of wall clock time to fully execute.

```console
carpdm.py -i clinically_relevant_amr.fna -p 20000 -t 16 -f /path/to/filter/fasta.fna
```

## CARPDM Output

Running the carpdm.py script as above will result in the following output structure, sorted here by creation order. Note that running with the --clean flag will only output the final probeset, oligo pool, amplification primers, and analysis plots.

```console
probe_design/
├── probes_input_no_comp.fna
|     Input fasta with complementary sequences between targets removed     
├── probes_basic_filter.fna
|     Naive tiled probes. Only unique probes with no perfect complements that
|     satisfy Tm requirements and lack ambiguous bases are included. If an
|     o-pool is being constructed, probes with an LguI cut site are also
|     removed.
├── probes_id_blast.txt
|     BLASTN results for basic filter probes against negative id filter
|     fasta/db (if supplied). 
├── probes_id_filter.fna
|     Basic filter probes that have <(probe_length * 0.625) identities against
|     anything in the id filter fasta/db (if supplied).
├── probes_self_blast.txt
|     BLASTN results for the basic filter or, if present, id filter probes
|     against themselves.
├── probes_self_filter.fna
|     Basic filter or, if present, id filter probes with redundant members
|     removed. The degree of redundancy maintained in the set is determined by
|     the probe number cutoff specified when the program is called. Included
|     with --clean.
├── probes_o_pool_amp_primers.fna
|     PCR amplification primers for the oligo pool. Included with --clean.
├── probes_o_pool_oligos.fna
|     Final oligo pool sequences, to be ordered if synthesizing probes
|     in-house. Included with --clean.
├── probes_final_blast.xml
|     BLASTN results of the probes against the input sequences, used to
|     construct summary data and plots.
├── probes_target_info.csv
|     Summary statistics of target and associated coverage by probes in set.
├── probes_target_probe_pairs.csv
|     Target:Probe pairs, each line is a probe that aligns to that target.
├── probes_probe_info.csv
|     Summary statistics of probes in set.
├── probes_count_info.csv
|     Number of probes remaining after each filtering step.
├── probes_max_id.txt
|     Maximum number of identites remaining between probes after the
|     redundancy filter.
└── probes_plots ( All included with --clean)
    ├── individual_target_coverages
    |     Directory containing coverage plots of all sequences used as input.
    ├── probe_gc.png
    ├── probe_gc.svg
    |     Violin plot showing probe GC content distribution.
    ├── probe_tm.png
    ├── probe_tm.svg
    |     Violin plot showing probe melt temp distribution.
    ├── probe_num_targets.png
    ├── probe_num_targets.svg
    |     Violin plot showing the distribution of the number of input targets
    |     per probe.
    ├── target_len.png
    ├── target_len.svg
    |     Violin plot showing target length distribution.
    ├── target_gc.png
    ├── target_gc.svg
    |     Violin plot showing target GC content distribution.
    ├── target_coverage_prop.png
    ├── target_coverage_prop.svg
    |     Violin plot showing target coverage proportion distribution.
    ├── target_probe_count.png
    ├── target_probe_count.svg
    |     Violin plot showing distribution of the number of probes per target.
    ├── target_coverage_depth.png
    ├── target_coverage_depth.svg
    |     Violin plot showing distribution of the target coverage depth.
    ├── target_coverage_stdev.png
    └── target_coverage_stdev.svg
          Violin plot showing distribution of the target coverage standard
          deviation.
```