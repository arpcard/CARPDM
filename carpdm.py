#!/usr/bin/env python3
from itertools import combinations_with_replacement, permutations
from Bio.SeqUtils import MeltingTemp as mt
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from Bio import SeqIO, SeqUtils
from Bio.Blast import NCBIXML
from hashlib import sha256
from Bio.Seq import Seq
from io import StringIO
from tqdm import tqdm
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
import argparse
import primer3
import shlex
import glob
import os
import re


def check_for_tax_data(
        blast_db: str
) -> bool:
    '''Check to make sure tax data is present if the nt database is the one
    being blasted against'''
    tax_dbs = ['taxdb.btd', 'taxdb.bti']
    blast_dir_contents = set(os.listdir(os.path.dirname(blast_db)))
    
    return all(tax_db in blast_dir_contents for tax_db in tax_dbs)


def count_seqs(
        input_fasta: str
) -> int:
    '''Count the number of sequences in a fasta file'''
    count = 0
    for record in SeqIO.parse(input_fasta, 'fasta'):
        count += 1

    return count


def collapse_slices(
        slice_lst: list[tuple[int, int]],
        result_lst: list = None
) -> list[tuple[int, int]]:
    '''Recursion that collapses a list of slices to a smaller list based on
    overlap betwen the tuples'''
    # Instantiate result lst. Used to pass results up the recursion/
    if result_lst is None:
        result_lst = []
    # Base case, where only one item is left in the list.
    if len(slice_lst) == 1:
        result_lst.append(slice_lst[0])
        return result_lst
    # Base case, where two items are left, but they get combined
    elif slice_lst[0][-1] >= slice_lst[1][0] and len(slice_lst) == 2:
        result_lst.append((slice_lst[0][0], max(slice_lst[0][1],
                                                slice_lst[1][1])))
        return result_lst
    # Recursive case, where the first two items in the list are combined
    # and there is still more to process
    elif slice_lst[0][-1] >= slice_lst[1][0]:
        result_lst.append((slice_lst[0][0], max(slice_lst[0][1],
                                                slice_lst[1][1])))
        return collapse_slices(slice_lst[2:], result_lst)
    # Recursive case, where the first two items cannot be combined
    # and there is still more to process
    else:
        result_lst.append(slice_lst[0])
        return collapse_slices(slice_lst[1:], result_lst)


def chop_record(
        input_record: SeqRecord,
        slices_to_rc: list[tuple[int, int]],
        probe_length: int
) -> list[SeqRecord]:
    '''Returns a list of SeqRecord objects that have been cut to remove
    complementary sequences specified by a list of slice values'''
    # First, convert slices into a more useable format, where they specify the
    # regions to be kept, instead of removed
    # Subtract 1 to account for zero-indexing used in python
    all_slices = [position for s in slices_to_rc for position in s]
    all_slices.insert(0, 0)
    all_slices.append(len(input_record.seq) + 1)
    # Instantiate slice to keep with 10bp of the complementary region,
    slices_to_keep = [(max(0, all_slices[i]), all_slices[i + 1])
                      for i in range(0, len(all_slices), 2)]
    # Instantiate list
    records = []
    # Cut out each slice. If it's a long enough sequence, add it to the record
    # list
    for s in slices_to_keep:
        working_record = input_record[s[0]:s[1]]
        if len(working_record.seq) >= probe_length:
            working_record.id = (working_record.id +f'_[{s[0]}:{s[1]}]')
            records.append(working_record)
    for s in slices_to_rc:
        working_record = input_record[s[0]:s[1]]
        if len(working_record.seq) >= probe_length:
            working_record.id = (working_record.id +f'_[{s[0]}:{s[1]}]')
            working_record.seq = working_record.seq.reverse_complement()
            records.append(working_record)
    
    return records


def cut_seqs(
        input_fasta: list[SeqRecord],
        slice_dict: dict[str, set[tuple[int, int]]],
        probe_length: int
) -> None:
    '''Return a list of SeqRecord objects with complementary regions removed'''
    # Instantiate set to quickly check membership
    slice_seqs = slice_dict.keys()
    new_fasta = []
    for record in input_fasta:
        # If the record needs to be sliced
        if record.id in slice_seqs:
            # Instantiate slices
            slices = slice_dict[record.id]
            # If there's more than one, istantiate a sorted list and check if
            # it can be collapsed
            old_slice_lst = sorted(list(slices), key = lambda x: (x[0], x[1]))
            if len(old_slice_lst) > 1:
                # Continue collapsing until no more progress is made
                while True:
                    new_slice_lst = collapse_slices(old_slice_lst)
                    if new_slice_lst == old_slice_lst or len(new_slice_lst
                                                             ) == 1:
                        final_slice_lst = new_slice_lst
                        break
                    old_slice_lst = new_slice_lst
            # If there's only one, proceed with that one
            else:
                final_slice_lst = old_slice_lst
            # Use the new slice lst to make new records that don't contain
            # those regions
            record_lst = chop_record(record, final_slice_lst, probe_length)
            new_fasta.extend(record_lst)
        # If it doesn't need to be sliced, append it for the next iteration
        else:
            new_fasta.append(record)
    # Return
    return new_fasta


def remove_complementary_targets(
        input_fasta: list[SeqRecord],
        output_dir: str,
        prefix: str,
        header: list[str],
        num_threads: int,
        probe_length: int
) -> None:
    '''Blasts a file against itself and removes any regions of complementarity.
    Returns a list of sorted Seq records'''
    # While loop to keep running until there are no more complementary hits
    while True:    
        # Make temporary BLAST db
        str_fasta = ''.join([record.format('fasta') for record in input_fasta])
        database = f'{output_dir}/self_db_temp'
        command = shlex.split(
            f'makeblastdb -out {database} -dbtype nucl -title self_db_temp')
        subprocess.run(command, input = str_fasta, text = True)
        # Perform BLAST
        outfmt = f'6 {" ".join(header)}'
        command = shlex.split(
            f'blastn -outfmt "{outfmt}" -num_threads {num_threads} '
            f'-db {database} -strand minus')
        blast_output = subprocess.run(command, text = True, input = str_fasta,
                                      capture_output = True).stdout
        self_blast_df = pd.read_csv(StringIO(blast_output), sep = '\t',
                                    names = header)
        # Pull out the sequence id with the most hits
        working_df = self_blast_df.loc[self_blast_df['nident'] >= 30]
        # If there are no more complementary hits, break the cycle and write
        # the fasta for later input
        if working_df.empty:
            SeqIO.write(input_fasta, f'{prefix}_input_no_comp.fna', 'fasta')
            break

        top_hit = working_df['qseqid'].mode().squeeze()
        if type(top_hit) == pd.Series:
            top_hit = top_hit.sort_values(ascending = False)[0]
        top_hit_df = working_df.loc[working_df['qseqid'] == top_hit]

        slice_dict = defaultdict(set)
        for t in top_hit_df[['sseqid', 'sstart', 'send']].itertuples():
            slice_dict[t[1]].add((min(t[2], t[3]), max(t[2], t[3])))
        input_fasta = cut_seqs(
            input_fasta = input_fasta,
            slice_dict = slice_dict,
            probe_length = probe_length
        )

    # Remove temporary database files
    db_files = [f'{output_dir}/{file}' for file in os.listdir(output_dir)
                if file.startswith('self_db_temp')]
    for file in db_files:
        os.remove(file)


def make_probes(
    design_fasta: str,
    probe_len: int,
    step: int
) -> set:
    '''Perform naive tiling across input fasta to generate a set of unique
    80bp probes'''
    all_probes = set()
    for record in SeqIO.parse(design_fasta, 'fasta'):
        # Instantiate uppercase string representation of the sequence
        sequence = record.seq.upper()
        # For every start position separated by step length up to the length
        # of the target - probe_length (doesn't add probes < probe_len) + 1
        # (accounts for exclusive terminal slice) add the probeseq starting at
        # that position to the probe list.
        all_probes.update({sequence[index:index + probe_len] for index in
                           range(0, len(record.seq) - probe_len + 1, step)})
        # If the step results in an uncovered portion at the end, append one
        # more probe that covers the last few bases
        if len(sequence) % step != 0:
            all_probes.add(sequence[-80:])
    
    return all_probes


def complementarity_filter(
        probe_set: set
) -> set:
    '''Remove probes with perfect complements in the probe_set'''
    # Make a sorted list for reproducibility
    probe_lst = sorted(list(probe_set))
    # Remove the probe sequence if its reverse complement is in the probe
    # set. By doing it this way, its reverse complement shouldn't be
    # removed, since it won't have a matching partner in the set.
    for probe in probe_lst:
        if probe.reverse_complement() in probe_set:
            probe_set.remove(probe)

    return probe_set
    

def basic_filter(
    probe_set: set,
    tm_low: int
) -> set:
    '''Filter probeset based on having a melting temperature >= tm_low'''
    prog = re.compile(r'[^ATGC]')
    # Filter out ambiguous bases and probes with too low a Tm
    filtered_probes = {
        probe for probe in probe_set if not prog.search(str(probe))
        and mt.Tm_NN(probe, nn_table = mt.R_DNA_NN1) > tm_low}
    
    return filtered_probes


def LguI_filter(probe_set: set) -> set:
    '''Remove probes with LguI cut sites; otherwise these will be cut during
    probe synthesis'''
    lgui_probes = {probe for probe in probe_set
                   if 'GAAGAGC' not in probe and 'GCTCTTC' not in probe}
    
    return lgui_probes


def write_to_file(
        probe_set: set,
        output: str
    ) -> None:
    '''Write to file, sorted based on sequence for reproducibility'''
    record_lst = []
    probe_lst = sorted(list(probe_set))
    for probe in probe_lst:
        record = SeqRecord(
            probe,
            id = f'probe_{sha256(str(probe).encode("utf-8")).hexdigest()}',
            name = '',
            description = ''
        )
        record_lst.append(record)
    SeqIO.write(record_lst, output, 'fasta')


def make_blast_db(
        filter_fasta: str,
        output: str
) -> None:
    # Make temporary BLAST db
    command = shlex.split(
        f'makeblastdb -in {filter_fasta} -out {output} -dbtype nucl')
    subprocess.run(command)


def blast_id(
        input_fasta: str,
        blast_db: str,
        header: list,
        num_threads: int,
        output: str
    ) -> None:
    '''Blasts probe sequences from a fasta against filter db then 
    saves to txt file'''
    os.environ['BLASTDB'] = os.path.dirname(blast_db)
    db_name = os.path.basename(blast_db)
    outfmt = f'6 {" ".join(header)}'
    command = shlex.split(f'blastn -query {input_fasta} -out {output} '
                          f'-outfmt "{outfmt}" -num_threads {num_threads} '
                          f'-db {db_name}')

    subprocess.run(command)


def filter_id_blast(
        id_blast_df: pd.DataFrame,
        probe_len: int
) -> set:
    '''Return set of probes that have disqualifying alignments
    from the id blast'''
    # Instantiate identity number cutoffs
    cutoff = probe_len * 0.625
    # Find the probes with > than that many identities
    probes_to_drop = set(id_blast_df.loc[(id_blast_df['nident'] < cutoff),
                                         'qseqid'])
    
    return probes_to_drop


def filter_nt_blast(
        id_blast_df: pd.DataFrame,
        probe_len: int
) -> set:
    '''Return set of probes that have disqualifying alignments
    from the id blast'''
    # Removed bacterial filter; if it has >80%id over <50nt,
    # it's not going to stick
    # Instantiate identity number cutoffs
    non_bac_high_cutoff = probe_len * 0.975
    non_bac_low_cutoff = probe_len * 0.625
    # Instantiate set of organisms for filter in the next step
    filter_set = set(['Archaea', 'Eukaryota', 'Viruses'])
    # Find the probes with eukaryotic, archaeal, or viral hits with number of
    # identities inside the cutoff range. By not removing those with 
    # nident > high cutoff, we avoid the problem of losing probes to erroneous
    # inclusion of a resistance target in one of these accessions.
    probes_to_drop = set(id_blast_df.loc[
        (id_blast_df['sskingdom'].isin(filter_set))
         & (id_blast_df['nident'] < non_bac_high_cutoff)
         & (id_blast_df['nident'] > non_bac_low_cutoff),
        'qseqid'])
    
    return probes_to_drop


def filter_out_probes(
        probe_fasta: str,
        probes_to_drop: set,
        output: str
) -> None:
    '''Remove probes listed in the set of probe names'''
    # Instantiate final list
    filtered_lst = []
    count = 0
    # For each record, if it's id isn't in the probe set, add it to the list
    for record in SeqIO.parse(probe_fasta, 'fasta'):
        count +=1
        if record.id not in probes_to_drop:
            filtered_lst.append(record)

    # Write resulting list to fasta file
    SeqIO.write(filtered_lst, output, 'fasta')


def blast_self(
        input_fasta: str,
        header: list,
        output: str,
        num_threads: int
    ) -> None:
    '''Constructs Blast DB from an input fasta, then blasts that fasta
    against the constructed db'''
    # Make temporary BLAST db
    output_dir = os.path.dirname(output)
    database = f'{output_dir}/self_db_temp'
    command = shlex.split(
        f'makeblastdb -in {input_fasta} -out {database} -dbtype nucl')
    subprocess.run(command)
    # Perform BLAST
    outfmt = f'6 {" ".join(header)}'
    command = shlex.split(
        f'blastn -query {input_fasta} -out {output} -outfmt "{outfmt}" '
        f'-num_threads {num_threads} -db {database}')
    subprocess.run(command)
    # Remove temporary database files
    db_files = [f'{output_dir}/{file}' for file in os.listdir(output_dir)
                if file.startswith('self_db_temp')]
    for file in db_files:
        os.remove(file)


def remove_complementary_probes(input_df: pd.DataFrame) -> set:
    '''Collapses probes with any complementarity to a single representative'''
    # Remove query:hit pairs that are for the same probe, and find entries
    # whose hit strans was 'Minus'. This indicates complementarity, not
    # identity
    rev_comp_df = input_df.loc[(input_df['qseqid'] != input_df['sseqid'])
                               # 30 cutoff based on melting temperature of RNA
                               # duplexes, and the assumption that there will
                               # often be gaps or substitutions that affect
                               # the melting temperature 
                               & (input_df['nident'] >= 30)
                               & (input_df['sstrand'] == 'minus')]
    # If any, start the process
    if rev_comp_df.empty == False:
        # Find the probes that had no queries or no hits,
        # add them to the list to remove
        probes_to_drop = set()
        query_set = set(rev_comp_df['qseqid'].to_list())    
        hit_set = set(rev_comp_df['sseqid'].to_list())
        probes_to_drop.update(hit_set ^ query_set)
        # Reinstantiate df without the probes with no queries/hits
        rev_comp_df = rev_comp_df.loc[
            (~rev_comp_df['qseqid'].isin(probes_to_drop))
            & (~rev_comp_df['sseqid'].isin(probes_to_drop))]
        # Get the number of hits to parse which probes to keep, and the
        # unique queries to iterate through
        hit_counts = rev_comp_df['qseqid'].value_counts()
        # This makes a list of queries sorted by the number of hits,
        # then by the query name to break ties. This allows the speed of a
        # greedy algorithm after the sort function, without needing to
        # iterate through every hit and check its value.
        queries = [k for k, v in sorted(hit_counts.items(),
                                        key = lambda item: (item[1],
                                                            item[0]))]
        # For each query
        for query in tqdm(queries):
            # If it hasn't already been slated to be dropped
            if query not in probes_to_drop:
                # Add the hits to the removal set
                hits = list(rev_comp_df.loc[rev_comp_df['qseqid'] == query,
                                            'sseqid'].unique())
                probes_to_drop.update(hits)
        return probes_to_drop
    else:
        print('No complementarity found')
        return set()


def get_num_sequences(input_fasta: str) -> int:
    '''Return the number of '>' in a fasta file'''
    with open(input_fasta) as infile:
        fasta = infile.read()
        return fasta.count('>')


def redundancy_filter(
        input_df: pd.DataFrame,
        cutoff: int,
        last_step_fasta: str,
        probes_to_drop: set = None
        ) -> set:
    '''Iteratively collapse probes with a certain number of identical
    residues to single representative. Check each time whether it's below
    a cutoff, then move to the next iteration if not'''
    # Get the number of current probes, and instantiate a working dataframe
    # that contains no probes aligned to themselves
    if probes_to_drop is None:
        probes_to_drop = set()
    nt_filter_probes = get_num_sequences(last_step_fasta)
    blast_df = input_df.loc[(input_df['qseqid'] != input_df['sseqid'])
                            &(input_df['sstrand'] == 'plus')]
    # Instantiate cutoff loop
    for id_cutoff in reversed(range(80)):
        # If you have less probes than the cutoff, break the iteration
        remaining_probes = nt_filter_probes - len(probes_to_drop)
        if remaining_probes < cutoff:
            break
        # Print progress
        print(f'Beginning redundancy filter with {id_cutoff} identities and '
              f'{remaining_probes} probes remaining')
        # Instantiate working df
        working_df = blast_df.loc[
            (blast_df['nident'] == id_cutoff)
            & (~blast_df['qseqid'].isin(probes_to_drop))
            & (~blast_df['sseqid'].isin(probes_to_drop))]
        # Find the probes that had no queries or no hits (orphans)
        hit_set = set(working_df['sseqid'].to_list())
        query_set = set(working_df['qseqid'].to_list())
        probes_to_drop.update(hit_set ^ query_set)
        # Reinstantiate working_df with no orphan probes, they're problematic
        working_df = working_df.loc[
            (~working_df['qseqid'].isin(probes_to_drop))
            & (~working_df['sseqid'].isin(probes_to_drop))]
        # Get the hit_counts
        hit_counts = working_df['qseqid'].value_counts()
        # This makes a list of queries sorted by the number of hits,
        # then by the query name to break ties. This allows the speed of a
        # greedy algorithm after the sort function, without needing to
        # iterate through every query and check its number of hits.
        queries = [k for k, v in sorted(hit_counts.items(),
                                        reverse = True,
                                        key = lambda item: (item[1],
                                                            item[0]))]
        for query in tqdm(queries):
            # If it hasn't already been slated to be dropped
            if query not in probes_to_drop:
                # Find the unique hits for that probe and append the query
                # sequence. This instantiates a list of the relevant probe
                # names to be compared.
                hits = list(working_df.loc[working_df['qseqid'] == query,
                                           'sseqid'].unique())
                probes_to_drop.update(hits)
    return probes_to_drop

### The following are scripts used if the set is being made for an o_pool

# T7_promoter:
#           5' - TAATACGACTCACTATAGGG - 3'
#       Transcription starts here ^
# Transcripts will contain three 5' Gs and the exact donwnstream
# sequence of this strand - opposite strand is the complement

# T7_promoter_complement:
#           3' - ATTATGCTGAGTGATATCCC - 5'
#       Transcription starts here ^
# Everything 5' of this will be used as template
# This is the complementary strand

# BspQI_cut_site:
#           5' - G C T C T T C N/N N N  - 3'
#           3' - C G A G A A G N N N N/ - 5'

# Transcription PCR molecule: 
# 5' - GCTAATACGACTCACTATAGGG|#probe_sequence#| NGAAGAGC|rev_comp_primer_sequence - 3'
# 3' - CGATTATGCTGAGTGATATCCC|probe_complement|/NCTTCTCG|final_primer_sequence### - 5'
# Other cut site is internal to the probe sequence, doesn't matter for
# positive strand, as it isn't used as template. This way, the probe
# will be properly sized + 3Gs at the 5' end

def make_primers(
        first_nts: str,
        last_nts: str,
        primer_length: int,
        primer2: str,
        output: str) -> None:
    '''Returns a sorted list of possible primers from all possible
    options given the constraints of the first and last nucleotides'''
    # Instantiate data types
    primer2_tm = primer3.calc_tm(primer2)
    count = 0
    record_list = []
    order_set = set()
    nts = 'ATGC'
    len_remaining = primer_length - len(first_nts) - len(last_nts)
    # Return all combinations of ATGC
    seqs = combinations_with_replacement(nts, len_remaining)
    for seq in seqs:
        # Return all unique permutations of all combinations. This is
        # the possible sequence space for these primers
        orders = set([''.join(order) for order in permutations(seq)])
        order_set.update(orders)
    # Sort the list
    sorted_orders = sorted(list(order_set))

    # Determine physical characteristics of the primers
    # Drop any that don't fit the bill
    for seq in sorted_orders:
        # Remove homopolymer repeats and LguI cut sites
        primer = first_nts + seq + last_nts
        if any([
            'AAA' in primer,
            'TTT' in primer,
            'GGG' in primer,
            'CCC' in primer,
            'GAAGAGC' in primer,
            'GCTCTTC' in primer[:-2]
        ]):
            continue
        # Calculate physical characteristics and compatibility with
        # primer 1
        primer_tm = primer3.calc_tm(primer)
        primer_hairpin = primer3.calc_hairpin(primer)
        primer_homodimer = primer3.calc_homodimer(primer)
        primer_heterodimer = primer3.calc_heterodimer(primer, primer2)
        if all([
            (primer2_tm - 0.1) < primer_tm < (primer2_tm + 0.1),
            primer_hairpin.tm < 20,
            primer_homodimer.tm < 20,
            primer_heterodimer.tm < 20
        ]):
            # Append primers that fit the physical characteristics
            count += 1
            record_list.append(
                SeqRecord(
                    Seq(primer),
                    id = (f'primer_'
                          f'{sha256(primer.encode("utf-8")).hexdigest()}'),
                    description = ''
                )
            )
    print(f'Finished contructing candidate primers, {count} found')
    # Write to file
    SeqIO.write(record_list, output, 'fasta')


def attach_seq(
        input_fasta: str,
        output: str,
        leading_sequence: str = '',
        lagging_sequence: str = ''
        ) -> None:
    '''Attach a sequence to either end of all sequences in a fasta
    file'''
    # Insatantiate variables
    record_list = []
    records = SeqIO.parse(input_fasta, 'fasta')
    # Add leading and lagging sequences
    for record in records:
        record.seq  = leading_sequence + record.seq + lagging_sequence
        record_list.append(record)
    # Rewrite as file
    SeqIO.write(record_list, output, 'fasta')


def blast_primers(
        query_file: str,
        sbjct_file: str,
        output: str,
        header: list
        ) -> None:
    '''Blasts primers against probes'''
    # BLAST file against the DB using settings for more sensitive
    # detection of short sequences
    outfmt = f'6 {" ".join(header)}'
    command = shlex.split(
        f'blastn -task "blastn-short" -query {query_file} -out {output} '
        f'-outfmt "{outfmt}" -subject {sbjct_file} -word_size 4 -evalue 50'
    )
    subprocess.run(command)


def parse_primer_blast(input_file: str, header: list):
    '''Returns the probe with the fewest identities to the probe_set'''
    # Instantiate variables
    primer_blast_df = pd.read_csv(input_file, sep = '\t', names = header)
    id_dict = {}
    # Determine the total number of identities for each primer
    for primer in primer_blast_df['qseqid'].unique():
        id_dict[primer] = primer_blast_df.loc[
            primer_blast_df['qseqid'] == primer, 'nident'
        ].sum()
    # Find the primer with the fewest identities to the probe set 
    final_primer = min(id_dict, key = id_dict.get)
    print(f'Finished parsing BLAST, {final_primer} performed best with '
          f'{id_dict[final_primer]} total identities '
          'to the original probeset')
    return final_primer


def fetch_sequence(name, input_fasta):
    '''Returns a sequence from a fasta file given the name'''
    records = SeqIO.parse(input_fasta, 'fasta')
    for record in records:
        if record.id == name:
            return str(record.seq)


def reverse_comp_sequence(seq):
    '''Returns the reverse complement of a sequence'''
    fwd = Seq(seq)
    rev_comp = str(fwd.reverse_complement())

    return rev_comp


def remove_temp(target_dir):
    to_remove = [f'{target_dir}/{file}' for file in os.listdir(target_dir)
                 if 'temp' in file]
    for file in to_remove:
        os.remove(file)

### These scripts are used for the analysis of the
## probes against the original fasta.

def blast(
        input_fasta: str,
        design_fasta: str,
        output: str
    ) -> None:
    '''Uses BLAST to compare the probe fasta to the sequences
    against which it was designed. Saves as xml to access match object.'''
    if os.path.exists(output) == False:
        command = shlex.split(
            f'blastn -query {input_fasta} -subject {design_fasta} '
            f'-out {output} -outfmt 5'
        )
        subprocess.run(command)


def make_dicts(
        probe_fasta: str,
        design_fasta: str
    ):
    '''Generate the necessary dictionaries for downstream processing'''
    # Workhorse; records the coverage of each
    # id in a target using a numpy array.
    coverage_dict = {}
    # Records the GC content of each target
    gc_dict = {}
    # Records the sequence of each target
    target_seq_dict = {}
    # Counts the number of target targets per probe
    probe_count_dict = {}
    # Counts the number of times a target is targeted by a probe
    target_count_dict = {}
    # Get target names
    target_records = SeqIO.parse(design_fasta, 'fasta')
    
    # Records the length and encodes an empty NumPy array of that
    # length for each target.
    # Also instantiates a zero in another dict to count the number of
    # times a target is targeted by a probe.
    for record in target_records:
        name = record.id.split('|')[-1]
        target_gc = SeqUtils.gc_fraction(record.seq)
        target_length = len(record.seq)
        coverage_dict[name] = np.zeros(target_length)
        gc_dict[name] = target_gc
        target_count_dict[name] = 0
        target_seq_dict[name] = str(record.seq)
    
    # Get probe names and instantiate a zero to count its number of
    # targets
    probe_records = SeqIO.parse(probe_fasta, 'fasta')
    for record in probe_records:
        name = record.id
        probe_count_dict[name] = 0

    # Return gc_dict and target_seq_dict separately from list, they
    # aren't needed immediately
    return (coverage_dict, probe_count_dict, target_count_dict,
            gc_dict, target_seq_dict)


def parse_blast(
    xml_file: str,
    coverage_dict: dict,
    probe_count_dict: dict,
    target_count_dict: dict,
    identity_cutoff: int,
    prefix: str
):
    '''Parse the blast xml file to determine coverage of each target,
    as well as which probes target which targets'''
    # Instantiate list to keep track of which probes target which targets
    pair_list = []
    with open(xml_file) as infile:
        # Parse every HSP individually, however only if they have more
        # identities than a give cutoff
        for record in NCBIXML.parse(infile):
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    # If there aren't enough identities,

                    # or there's a gap >1nt, don't consider
                    if hsp.identities > identity_cutoff and not any([
                        re.search(r'--+', hsp.sbjct) != None,
                        re.search(r'--+', hsp.query) != None
                        ]):
                        # Extract target name from alignment title
                        target = alignment.title.split('|')[-1].split()[0]
                        # Each hit means another probe that targets
                        target_count_dict[target] += 1
                        # Extract probe name from query
                        probe = record.query
                        # Each hit means another targeted target
                        probe_count_dict[probe] += 1
                        # Append the target:probe pair to list
                        pair_list.append((target, probe))
                        # target length = NumPy array size
                        target_len = coverage_dict[target].size
                        # Convert match characters ('|' for match,
                        # ' ' for no match) to string of space-
                        # delimited integers ('1' for match,
                        # '0' for no match) 
                        matches = ' '.join(
                            hsp.match.replace('|', '1').replace(' ', '0')
                        )
                        # Convert integer string to NumPy array
                        match_array = np.fromstring(
                            matches, dtype=int, sep = ' '
                        )
                        # If there's a gap in the subject, remove
                        # the corresponding base from the match
                        # array. This fixes a problem that arises
                        # with improper padding.
                        if '-' in hsp.sbjct:
                            sbjct_gaps = ' '.join(
                                re.sub(
                                    r'[ATGC]',
                                    '0',
                                    hsp.sbjct.replace('-', '1')
                                )
                            )
                            sbjct_gap_array = np.fromstring(
                                sbjct_gaps, dtype=int, sep = ' '
                            )
                            sbjct_gap_position = np.where(
                                sbjct_gap_array == 1
                            )
                            match_array = np.delete(
                                match_array, sbjct_gap_position
                            )
                        # Use the preceding information to pad the
                        # match array with sufficient zeros to make
                        # it the same size as the corresponding
                        # target array.
                        hsp_length = match_array.size
                        # If hsp strands aren't equivalent, indicates the
                        # start position is on the opposite side.
                        if hsp.strand[0] == hsp.strand[1]:
                            start = hsp.sbjct_start
                        else:
                            start = hsp.sbjct_end
                        # Pad array with zeros
                        pad_left = start - 1
                        pad_right = target_len - (pad_left + hsp_length)
                        padded_array = np.pad(
                            match_array,
                            (pad_left, pad_right),
                            mode = 'constant',
                            constant_values = (0,0)
                        )
                        # Add the padded array to the corresponding
                        # target coverage array, such that each
                        # nucleotide position with a match will
                        # increment by 1. This represents the 
                        # coverage at each position.
                        coverage_dict[target] += padded_array
    # Write probe:target pairs to a summary csv
    pair_df = pd.DataFrame(pair_list, columns = ['target', 'probe'])
    pair_df.to_csv(f'{prefix}_target_probe_pairs.csv', index = False)
    return coverage_dict, target_count_dict, probe_count_dict


def parse_targets(
        coverage_dict: dict,
        target_count_dict: dict,
        target_gc_dict: dict,
        target_seq_dict: dict,
        prefix: str
    ):
    '''Extracts relevant target info and writes to csv'''
    # Instantiate summary df
    colnames = ['target_name', 'length', 'gc', 'sequence', 'probe_count',
                'proportion_covered', 'mean_depth','std_dev']
    target_info = pd.DataFrame(columns = colnames)
    # Parse relevant info
    for target in coverage_dict:
        probe_count = target_count_dict[target]
        array = coverage_dict[target]
        length = array.size
        gc = target_gc_dict[target]
        seq = target_seq_dict[target]
        mean_depth = np.mean(array)
        std_dev = np.std(array)
        num_covered = np.count_nonzero(array)
        num_positions = array.size
        proportion_covered = num_covered/num_positions
        # Instantiate list of relevant info
        info_list = [
            target, length, gc, seq, probe_count,
            proportion_covered, mean_depth, std_dev
        ]
        # Append list to summary df
        target_info.loc[len(target_info)] = info_list
    # Write summary df to csv
    target_info.to_csv(f'{prefix}_target_info.csv', index = False)
    
    
def parse_probes(
        probe_fasta: str,
        probe_count_dict: dict,
        prefix: str
    ) -> None:
    '''Extracts relevant probe info and writes to csv'''
    # Instantiate summary df
    colnames = ['probe_id', 'gc', 'tm', 'num_targets', 'sequence']
    probe_info = pd.DataFrame(columns = colnames)
    # Read in probe records
    probe_records = SeqIO.parse(probe_fasta, 'fasta')
    # For each probe, parse out relevant info
    for record in probe_records:
        probe_id = record.id
        num_targets = probe_count_dict[probe_id]
        gc = SeqUtils.gc_fraction(record.seq)
        tm = mt.Tm_NN(record.seq)
        sequence = str(record.seq)
        # Assemble info list, append to existing df.
        info_list = [ probe_id, gc, tm, num_targets, sequence]
        # Append list to summary df
        probe_info.loc[len(probe_info)] = info_list
    
    # Write to summary csv
    probe_info.to_csv(f'{prefix}_probe_info.csv', index = False)


def get_probe_names(input_fasta: str) -> set:
    '''Extract probe names from a fasta file, return them as a set'''
    with open(input_fasta) as infile:
        probe_names = set([
            record.description for record in SeqIO.parse(infile, 'fasta')
        ])
    return probe_names
    

def find_remaining_ID(
        probe_set: str,
        self_blast: str,
        header: list,
        prefix: str
) -> None:
    '''Determines the max remaining number of identities between probes
    in the final probeset'''
    #Find the max number of identities
    probe_names = get_probe_names(probe_set)
    # Read in the self_blast tsv
    self_blast_df = pd.read_csv(self_blast, sep = '\t', names = header)
    # Determine the max remaining number of identities in the probeset
    max_id = self_blast_df.loc[
        self_blast_df['qseqid'].isin(probe_names)
        & self_blast_df['sseqid'].isin(probe_names)
        & (self_blast_df['qseqid'] != self_blast_df['sseqid']),
        'nident'].max()
    print(f'Max identities between remaining probes in the set is {max_id}')
    # Save to .txt file
    with open(f'{prefix}_max_id.txt', 'w+') as outfile:
        outfile.write(str(max_id))


def count_probes(prefix: str) -> None:
    '''Records the number of probes after each filter'''
    # Find fastas
    fasta_list = glob.glob(f'{prefix}*.fna')
    # Instantiate summary df
    col_names = ['filter', 'count']
    count_df = pd.DataFrame(columns = col_names)
    
    # Count the number of probes in each fasta
    for fasta in fasta_list:
        count = 0
        records = SeqIO.parse(fasta, 'fasta')
        for _ in records:
            count += 1
        info_list = [fasta, count]
        count_df.loc[len(count_df)] = info_list

    # Write summary file
    count_df.to_csv(f'{prefix}_count_info.csv', index = False)


def plot_individual_coverages(
        coverage_dict: dict[str, np.ndarray],
        plot_output_dir: str
):
    for target in coverage_dict:
        cov_array = coverage_dict[target]
        fig, ax = plt.subplots()
        sns.lineplot(x = np.arange(len(cov_array)), y = cov_array, ax = ax)
        ax.set_title(f'{target} Coverage')
        plt.savefig(
            f'{plot_output_dir}/{target.replace("/", "|")}_coverage.png',
            dpi = 150)
        plt.close(fig)


def plot_probe_gc(
    probe_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot probe set gc content distribution as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = probe_df['gc'], ax = ax)
    ax.set_title('Probe Set GC Content Distribution')
    plt.savefig(f'{plot_output_dir}/probe_gc.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/probe_gc.svg')
    plt.close(fig)


def plot_probe_tm(
    probe_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot probe set melting temperature distribution as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = probe_df['tm'])
    ax.set_title('Probe Set Melt Temp Distribution')
    plt.savefig(f'{plot_output_dir}/probe_tm.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/probe_tm.svg')
    plt.close(fig)


def plot_probe_num_targets(
    probe_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot number of targets per probe as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = probe_df['num_targets'])
    ax.set_title('Probe Set Target Number Distribution')
    plt.savefig(f'{plot_output_dir}/probe_num_targets.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/probe_num_targets.svg')
    plt.close(fig)


def plot_target_len(
    target_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot distribution of target lengths as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = target_df['length'])
    ax.set_title('Target Length Distribution')
    plt.savefig(f'{plot_output_dir}/target_len.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/target_len.svg')
    plt.close(fig)


def plot_target_gc(
    target_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot distribution of target lengths as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = target_df['gc'])
    ax.set_title('Target GC Content Distribution')
    plt.savefig(f'{plot_output_dir}/target_gc.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/target_gc.svg')
    plt.close(fig)


def plot_target_probe_counts(
    target_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot distribution of target probe counts as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = target_df['probe_count'])
    ax.set_title('Target Probe Count Distribution')
    plt.savefig(f'{plot_output_dir}/target_probe_count.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/target_probe_count.svg')
    plt.close(fig)


def plot_target_coverage_prop(
    target_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot distribution of target lengths as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = target_df['proportion_covered'])
    ax.set_title('Target Coverage Proportion Distribution')
    plt.savefig(f'{plot_output_dir}/target_coverage_prop.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/target_coverage_prop.svg')
    plt.close(fig)


def plot_target_coverage_depth(
    target_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot distribution of target coverage depth as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = target_df['mean_depth'])
    ax.set_title('Target Coverage Depth Distribution')
    plt.savefig(f'{plot_output_dir}/target_coverage_depth.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/target_coverage_depth.svg')
    plt.close(fig)


def plot_target_coverage_stdev(
    target_df: pd.DataFrame,
    plot_output_dir: str
) -> None:
    '''Plot distribution of target coverage std deviation as violin plot'''
    fig, ax = plt.subplots()
    sns.violinplot(y = target_df['std_dev'])
    ax.set_title('Target Coverage Std Deviation Distribution')
    plt.savefig(f'{plot_output_dir}/target_coverage_stdev.png', dpi = 300)
    plt.savefig(f'{plot_output_dir}/target_coverage_stdev.svg')
    plt.close(fig)


def clean(
        output_dir: str
) -> None:
    '''Remove unnecessary files when --clean is specified'''
    all_files = [f'{output_dir}/{file}' for file in os.listdir(output_dir)]
    to_remove = [f for f in all_files if not f.endswith('self_filter.fna')
                 and not f.endswith('plots') and not 'o_pool' in f]

    for file in to_remove:
        os.remove(file)
    


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = (
            'Design probes against an input fasta according to design '
            'parameters, then pass through several filters to remove '
            'off-target enrichment (id filter) and probe redundancy '
            '(self filter)'
        )
    )
    parser.add_argument(
        '-i', '--input_fasta', required = True,
        help = 'Path to the fasta that contains sequences against which '
        'probes are to be designed. Required'
    )
    parser.add_argument(
        '-p', '--probe_num_cutoff', default = 42000, type = int,
        help = 'Maximum number of final probes allowed in the probeset. '
        'The redundancy filter will iterate until it reaches below this '
        'cutoff. Default = 42000'
    )
    parser.add_argument(
        '-s', '--tiling_step', default = 4, type = int,
        help = 'Number of bases by which to offset initial probe tiling. '
        'Higher values make computation less intense, though may result in '
        'sparser coverage. Default = 4'
    )
    parser.add_argument(
        '-l', '--probe_length', default = 80, type = int,
        help = 'Probe length to create. Note that if >80 (default), probes '
        'will be too long for the base price of a Twist oligo pool after '
        'addition of amplification and transcription primers for in-house '
        'synthesis. Default = 80'
    )
    parser.add_argument(
        '-m', '--melt_temp', default = 50, type = int,
        help = 'Minimum melting temperature for probes to be considered '
        'during basic filter. Default = 50'
    )
    parser.add_argument(
        '-b', '--basename',
        help = 'File basename. Default = probes',
        default = 'probes'
    )
    parser.add_argument(
        '-o', '--output_dir',
        help = 'Output directory. Default = probe_design',
        default = 'probe_design'
    )
    filter_group = parser.add_mutually_exclusive_group()
    filter_group.add_argument(
        '-f', '--filter_fasta',
        help = 'Fasta file to filter against. This creates a blast database '
        'from the provided file to filter probes against after the basic '
        'filter. Probes with over >(probe_length * 0.625) identities '
        'against anything in this fasta file will be removed' 
    )
    filter_group.add_argument(
        '-d', '--filter_db',
        help = 'Premade blast database to filter against. Probes with '
        '>(probe_length * 0.625) identities against anything in this '
        'fasta file will be removed, unless the specified database is the '
        'nt database (basename == "nt"). In that case it will only remove '
        'probes that align with that condition to non-bacterial accessions, '
        'that also do not have >(probe_length * 0.975) identities'
    )
    parser.add_argument(
        '-t', '--num_threads', type = int, default = 1,
        help = 'Number of threads to be used during BLAST searching'
    )
    parser.add_argument(
        '-n', '--no_o_pool', action = 'store_true',
        help = 'If included, omits steps to design template oligo pool. '
        'Only include if ordering RNA probes directly from the manufacturer.'
    )
    parser.add_argument(
        '-c', '--clean', action = 'store_true',
        help = 'If included, removes all files except for final baits and '
        'oligo pools, as well as analysis plots'
    )
    parser.add_argument(
        '-v', '--version', action = 'version', version = '%(prog)s v1.0.1'
    )

    return parser.parse_args()


def main():
    ############################### Parse Args ###############################
    args = parse_arguments()
    design_fasta = args.input_fasta
    probe_num_cutoff = args.probe_num_cutoff
    probe_len = args.probe_length
    basename = args.basename
    output_dir = args.output_dir
    step = args.tiling_step
    no_o_pool = args.no_o_pool
    num_threads = args.num_threads
    melt_temp = args.melt_temp
    filter_fasta = args.filter_fasta
    filter_db = args.filter_db
    clean_files = args.clean

    # Preliminary checks
    if filter_db:
        if os.path.basename(filter_db) == 'nt':
            assert check_for_tax_data(filter_db), 'No tax data in nt blast db directory'
            print('taxonomy data present in blast nt db directory')

    # Make output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Join output directory and file basename for prefix
    prefix = os.path.join(output_dir, basename)
    
    # Instantiate BLAST results header
    header = ('qseqid sseqid length pident nident staxid sstrand qstart qend '
              'sstart send qlen slen').split()
    
    ############################# Probe Design ###############################
    
    ## QC for input
    # Perform iterative blast of input queries to determine if there are
    # any complementary regions, and remove them if so
    if not os.path.exists(f'{prefix}_input_no_comp.fna'):
        input_fasta = list(SeqIO.parse(design_fasta, 'fasta'))
        remove_complementary_targets(
            input_fasta = input_fasta,
            output_dir = output_dir,
            prefix = prefix,
            header = header,
            num_threads = num_threads,
            probe_length = probe_len
        )
    else:
        print('Probes seem to have already undergone initial construction, '
              'skipping to basic filter')

    ## Basic Filters
    if not os.path.exists(f'{prefix}_basic_filter.fna'):
        print('Beginning probe construction')
        probes = make_probes(
            design_fasta = f'{prefix}_input_no_comp.fna',
            probe_len = probe_len,
            step = step
        )
        write_to_file(probes, f'{prefix}_no_filter.fna')
        print(f'{len(probes)} probes constructed during naive tiling')
        no_comp_probes = complementarity_filter(probes)
        print(f'{len(no_comp_probes)} probes remaining after removing '
              'perfect complements')
        filtered_probes = basic_filter(no_comp_probes, melt_temp)
        if not no_o_pool:
            filtered_probes = LguI_filter(filtered_probes)
        basic_probe_num = len(filtered_probes)
        print(f'{basic_probe_num} probes remaining after basic filter')
        write_to_file(filtered_probes, f'{prefix}_basic_filter.fna')
    else:
        print('Probes seem to have already undergone basic filter, '
              'skipping to BLAST id filter')
        basic_probe_num = count_seqs(f'{prefix}_basic_filter.fna')
        


    ## Blast ID Filter
    # If the id filter hasn't already been done
    if not os.path.exists(f'{prefix}_id_filter.fna'):
        # If there's a fasta to filter against, make a blast database
        if filter_fasta:
            make_blast_db(
                filter_fasta = filter_fasta,
                output = f'{prefix}_filter_db_temp'
            )
            # This overwrites the previous value of None if a filter fasta was
            # included, making it truthy.
            filter_db = f'{prefix}_filter_db_temp'
        if filter_db:
            # If blasting against nt database, only remove hits against
            # non-bacterial accessions
            if os.path.basename(filter_db) == 'nt':
                print('nt database detected for id filter, changing rules to '
                      'only exclude non-bacterial hits')
                # Reinstantiate header, start, end positions and pid no longer matter
                header = 'qseqid sseqid length nident staxid sstrand sskingdom'.split()
                blast_id(
                    input_fasta = f'{prefix}_basic_filter.fna',
                    blast_db = filter_db,
                    header = header,
                    num_threads = num_threads,
                    output = f'{prefix}_id_blast.txt'
                )
                id_blast_df = pd.read_csv(
                    f'{prefix}_id_blast.txt',
                    sep = '\t',
                    names = header
                )
                id_drop_probes = filter_nt_blast(
                    id_blast_df=id_blast_df,
                    probe_len=probe_len
                )
            # Otherwise, use standard rules
            else:
                header = 'qseqid sseqid length nident'.split()
                blast_id(
                    input_fasta = f'{prefix}_basic_filter.fna',
                    blast_db = filter_db,
                    header = header,
                    num_threads = num_threads,
                    output = f'{prefix}_id_blast.txt'
                )
                id_blast_df = pd.read_csv(
                    f'{prefix}_id_blast.txt',
                    sep = '\t',
                    names = header
                )
                id_drop_probes = filter_id_blast(
                    id_blast_df=id_blast_df,
                    probe_len=probe_len
                )
            filter_out_probes(
            probe_fasta = f'{prefix}_basic_filter.fna',
            probes_to_drop = id_drop_probes,
            output = f'{prefix}_id_filter.fna'
            )
            id_probe_num = basic_probe_num - len(id_drop_probes)
            print(f'{id_probe_num} probes remaining after id filter')
            print(f'{len(id_drop_probes)} probes removed during id filter')
            last_step_fasta = f'{prefix}_id_filter.fna'
        else:
            print('No id filter specified. Skipping to redundancy filter')
            last_step_fasta = f'{prefix}_basic_filter.fna'
    else:
        print('Probes seem to have already undergone BLAST id filter, '
              'skipping to redundancy/complementarity filter')
        last_step_fasta = f'{prefix}_id_filter.fna'
    id_probe_num = count_seqs(last_step_fasta)


    ## Complementarity/Redundancy Filter
    header = 'qseqid sseqid length nident sstrand'.split()
    if not os.path.exists(f'{prefix}_self_filter.fna'):
        blast_self(
            input_fasta = last_step_fasta,
            header = header,
            output = f'{prefix}_self_blast.txt',
            num_threads = num_threads
        )
        self_blast_df = pd.read_csv(
            f'{prefix}_self_blast.txt',
            sep = '\t',
            names = header
        )
        complementary_probes = remove_complementary_probes(self_blast_df)
        print(f'{len(complementary_probes)} probes removed during '
              'complementarity filter')
        redundant_probes = redundancy_filter(
            input_df = self_blast_df,
            cutoff = probe_num_cutoff,
            last_step_fasta = last_step_fasta,
            probes_to_drop = complementary_probes
        )
        filter_out_probes(
            probe_fasta = last_step_fasta,
            probes_to_drop = redundant_probes,
            output = f'{prefix}_self_filter.fna'
        )
        red_probe_num = id_probe_num - len(redundant_probes)
        print(f'{red_probe_num} probes remaining after redundancy filter')
        print(f'{len(redundant_probes)} probes removed during redundancy filter')
    else:
        print('Probes seem to have already undergone redundancy filter, '
              'skipping to amplification primer design and selection')
        red_probe_num = count_seqs(f'{prefix}_self_filter.fna')

    ############################ O-Pool Design #############################

    ## Oligo pool design ###
    # The following scripts are only run if for_o_pool == true, they add
    # T7 and another primer sequence to each end for amplification
    if not no_o_pool and not os.path.exists(
        f'{prefix}_o_pool_amp_primers.fna'):
        # Added a 5'GC to T7 primer to bump up the Tm a little
        T7_primer = 'GCTAATACGACTCACTATAGGG'
        # Add terminal GC clamp and BspQI to second primer.
        # Narrows down possibilities to a reasonable number
        # to minimize run time
        # Find all reasonable primer candidates
        print('Beginning construction of candidate primers')
        make_primers(
            first_nts = '',
            last_nts = 'GCTCTTCG',
            primer_length = 18,
            primer2 = T7_primer,
            output = f'{prefix}_possible_primers_temp.fasta'
        )
        # Attach the T7 primer sequence to the 5' end of the probes
        attach_seq(
            input_fasta = f'{prefix}_self_filter.fna',
            output = f'{prefix}_T7_probes_temp.fasta',
            leading_sequence = T7_primer
        )
        print('Finished constructing sequences with appended T7 primer')
        # Blast the primers against the probes with attached T7 primer
        blast_primers(
            query_file = f'{prefix}_possible_primers_temp.fasta',
            sbjct_file = f'{prefix}_T7_probes_temp.fasta',
            output = f'{prefix}_primer_blast_temp.txt',
            header = header
        )
        print('Finished BLASTing primers against probe set')
        # Find the primer with the fewest identities to the probes
        final_primer = parse_primer_blast(
            f'{prefix}_primer_blast_temp.txt',
            header = header
        )
        # Find the primer sequence from its name
        final_primer_sequence = fetch_sequence(
            name = final_primer, 
            input_fasta = f'{prefix}_possible_primers_temp.fasta'
        )
        # Reverse complement the primer
        primer_rev_comp = reverse_comp_sequence(final_primer_sequence)
        # Attach the reverse complement of the primer to the 3' end of the
        # probes with T7 primer attached, write to fasta
        attach_seq(
            input_fasta = f'{prefix}_T7_probes_temp.fasta',
            output = f'{prefix}_o_pool_oligos.fna',
            lagging_sequence = primer_rev_comp
        )
        print('Finished constructing full o-pool sequences')
        # Write primer sequences to fasta
        with open(f'{prefix}_o_pool_amp_primers.fna', 'w') as outfile:
            outfile.write(
                f'>T7_primer\n'
                f'{T7_primer}\n'
                f'>{final_primer}\n'
                f'{final_primer_sequence}\n'
            )
    else:
        print('Probes seem to already have amplification primers appended, '
              'skipping to analysis')

    ########################### BLAST Analysis ############################

    # Blast desired file
    if not os.path.exists(f'{prefix}_plots/target_coverage_stdev.svg'):
        print('Beginning blast analysis against input sequences')
        blast(
            input_fasta = f'{prefix}_self_filter.fna',
            design_fasta = design_fasta,
            output = f'{prefix}_final_blast.xml'
        )

        # Construct summary dicts
        (coverage_dict, probe_count_dict, target_count_dict,
        target_gc_dict, target_seq_dict) = make_dicts(
            probe_fasta = f'{prefix}_self_filter.fna',
            design_fasta = design_fasta
        )

        # Parse the blast output, also writes target_probe_pairs to csv
        print('Blast finished, beginning parsing')
        coverage_dict, target_count_dict, probe_count_dict = parse_blast(
            xml_file = f'{prefix}_final_blast.xml',
            coverage_dict = coverage_dict,
            probe_count_dict = probe_count_dict,
            target_count_dict = target_count_dict,
            identity_cutoff = 50,
            prefix = prefix
        )
        # Generate and write probe and target info to .csv files.
        print(
            'Parsing finished, beginning summarization and writing to .csv files')
        parse_targets(
            coverage_dict = coverage_dict,
            target_count_dict = target_count_dict,
            target_gc_dict = target_gc_dict,
            target_seq_dict = target_seq_dict,
            prefix = prefix
        )
        parse_probes(
            probe_fasta = f'{prefix}_self_filter.fna',
            probe_count_dict = probe_count_dict,
            prefix = prefix
        )
        # Find the max remaining identities between probes
        print('Summary csv files created, finding remaining identities')
        find_remaining_ID(
            probe_set = f'{prefix}_self_filter.fna',
            self_blast = f'{prefix}_self_blast.txt',
            header = header,
            prefix = prefix
        )
        # Finally, count the probes in each
        print('Remaining identities determined, counting probes in each set')
        count_probes(prefix)

        ## Plotting stats
        print('Beginning statistic plotting')
        # Make plot directories
        if not os.path.exists(f'{prefix}_plots'):
            os.mkdir(f'{prefix}_plots')
        if not os.path.exists(f'{prefix}_plots/individual_target_coverages'):
            os.mkdir(f'{prefix}_plots/individual_target_coverages')

        # Start with individual coverages
        print('Plotting individual target coverages')
        plot_individual_coverages(
            coverage_dict = coverage_dict,
            plot_output_dir = f'{prefix}_plots/individual_target_coverages'
        )

        # Read in probe df
        print('Plotting probe info')
        probe_df = pd.read_csv(f'{prefix}_probe_info.csv')
        plot_probe_gc(
            probe_df = probe_df,
            plot_output_dir = f'{prefix}_plots'
        )
        plot_probe_tm(
            probe_df = probe_df,
            plot_output_dir = f'{prefix}_plots'
        )
        plot_probe_num_targets(
            probe_df = probe_df,
            plot_output_dir = f'{prefix}_plots'
        )

        # Read in target_info_df
        print('Plotting target info')
        target_df = pd.read_csv(f'{prefix}_target_info.csv')
        plot_target_len(
            target_df = target_df,
            plot_output_dir = f'{prefix}_plots'
        )
        plot_target_gc(
            target_df = target_df,
            plot_output_dir = f'{prefix}_plots'
        )
        plot_target_probe_counts(
            target_df = target_df,
            plot_output_dir = f'{prefix}_plots'
        )
        plot_target_coverage_prop(
            target_df = target_df,
            plot_output_dir = f'{prefix}_plots'
        )
        plot_target_coverage_depth(
            target_df = target_df,
            plot_output_dir = f'{prefix}_plots'
        )
        plot_target_coverage_stdev(
            target_df = target_df,
            plot_output_dir = f'{prefix}_plots'
        )

        ## Clean up
        print('Finished plotting, removing temporary files')
        remove_temp(output_dir)
        if clean_files:
            clean(output_dir)
        print('Done!')
    else:
        print('Analysis seems to have already been completed')

if __name__ == '__main__':
    main()
