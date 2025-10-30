#!/usr/bin/env python3
"""
cdna_translate.py

Translate cDNA (coding DNA) to amino acid sequence.

Usage examples:
    python cdna_translate.py -s ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    python cdna_translate.py -i input.fasta -o protein.fasta --frame 0
    python cdna_translate.py -s ATG... --frame 1 --reverse

Features:
 - Accepts plain sequence via -s/--sequence or a FASTA file via -i/--input
 - Accepts both T and U
 - Frame can be 0, 1, or 2
 - Optionally translate reverse complement with --reverse
 - Outputs FASTA or plain sequence; wraps FASTA lines at 60 chars
"""

import argparse
import textwrap
import sys
from typing import Dict

CODON_TABLE: Dict[str, str] = {
    # Standard genetic code
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def parse_fasta(file_path: str) -> Dict[str, str]:
    """Simple FASTA parser returning dict header->sequence."""
    seqs = {}
    header = None
    seq_lines = []
    with open(file_path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            seqs[header] = ''.join(seq_lines)
    return seqs

def clean_seq(s: str) -> str:
    """Uppercase, remove whitespace and convert U->T"""
    return ''.join(s.split()).upper().replace('U', 'T')

def revcomp(seq: str) -> str:
    complement = str.maketrans('ATGCatgc', 'TACGtacg')
    return seq.translate(complement)[::-1]

def translate_seq(dna: str, frame: int = 0) -> str:
    """Translate DNA (T-based) starting at frame 0/1/2. Returns AA string with '*' for stops."""
    dna = dna.upper()
    aa = []
    for i in range(frame, len(dna) - 2, 3):
        codon = dna[i:i+3]
        if len(codon) < 3:
            break
        if any(base not in 'ATGC' for base in codon):
            # if ambiguous base present, use X
            aa.append('X')
            continue
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)

def main(argv=None):
    p = argparse.ArgumentParser(description="Translate cDNA to amino acid sequence.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--sequence', help='DNA sequence string (T or U allowed)')
    group.add_argument('-i', '--input', help='Input FASTA file with one or more DNA sequences')
    p.add_argument('-o', '--output', help='Output file (defaults to stdout)')
    p.add_argument('--frame', type=int, choices=(0,1,2), default=0, help='Reading frame (0,1,2). default=0')
    p.add_argument('--reverse', action='store_true', help='Translate reverse complement instead of given strand')
    p.add_argument('--fasta', action='store_true', help='Write output in FASTA format (header required for -s)')
    p.add_argument('--wrap', type=int, default=60, help='Line wrap width for FASTA output (default 60)')
    args = p.parse_args(argv)

    outputs = []

    if args.sequence:
        seq = clean_seq(args.sequence)
        if args.reverse:
            seq = revcomp(seq)
        prot = translate_seq(seq, frame=args.frame)
        header = 'sequence_translation'
        outputs.append((header, prot))
    else:
        # input FASTA: translate each record
        seqs = parse_fasta(args.input)
        if not seqs:
            print(f"No FASTA records found in {args.input}", file=sys.stderr)
            sys.exit(2)
        for hdr, s in seqs.items():
            s_clean = clean_seq(s)
            if args.reverse:
                s_clean = revcomp(s_clean)
            prot = translate_seq(s_clean, frame=args.frame)
            outputs.append((hdr, prot))

    # write output
    out_fh = open(args.output, 'w') if args.output else sys.stdout
    try:
        for hdr, prot in outputs:
            if args.fasta:
                out_fh.write(f'>{hdr}\n')
                out_fh.write('\n'.join(textwrap.wrap(prot, args.wrap)))
                out_fh.write('\n')
            else:
                # plain text: header + sequence on same line (tab separated)
                if len(outputs) == 1:
                    out_fh.write(prot + '\n')
                else:
                    out_fh.write(hdr + '\t' + prot + '\n')
    finally:
        if args.output:
            out_fh.close()

def find_longest_m_peptide(seq, translate_seq):
    """
    Finds the frame with the longest peptide starting with 'M' from a DNA sequence.
    
    Args:
        seq (str): DNA sequence.
        translate_seq (function): Function that translates a nucleotide sequence into amino acids.
    
    Returns:
        list: [best_frame, longest_length, stop_index]
              If no peptide starting with 'M' is found, returns [None, 0, None].
    """
    best_frame = None
    longest_length = 0
    best_stop_index = None

    for frame in range(3):
        # Translate the sequence starting from this frame
        aa_seq = translate_seq(seq[frame:], frame=0)
        aa_list = [frag for frag in aa_seq.split('*') if frag]

        # Keep only peptides starting with M
        m_peptides = [frag for frag in aa_list if frag.startswith('M')]

        if not m_peptides:
            continue

        # Find the longest M-starting peptide
        longest = max(m_peptides, key=len)
        length = len(longest)

        # Locate in the full amino acid sequence
        start_index = aa_seq.find(longest)
        stop_index = aa_seq.find('*', start_index)

        # Update if this is the longest so far
        if length > longest_length:
            best_frame = frame
            longest_length = length
            best_stop_index = stop_index if stop_index != -1 else None

    return [best_frame, longest_length, best_stop_index]


# parse_hgTables.py
def parse_hg_tables(input_file="hgTables.txt", output_file="hgTables_parsed.txt"):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        name = None
        sequence_parts = []

        for line in infile:
            line = line.strip()
            if not line:
                continue  # skip empty lines

            if line.startswith(">"):
                # Write previous sequence if we have one
                if name and sequence_parts:
                    sequence = "".join(sequence_parts)
                    outfile.write(f"{name}\t{sequence}\n")

                # Extract first token after '>'
                name = line[1:].split()[0]

                # Remove unwanted prefixes
                prefixes_to_remove = [
                    "mm39_ct_allisoforms_2790_",
                    "CACNA1C_"
                ]
                for prefix in prefixes_to_remove:
                    if name.startswith(prefix):
                        name = name[len(prefix):]

                sequence_parts = []  # reset sequence buffer
            else:
                # Append sequence lines
                sequence_parts.append(line)

        # Write the last sequence
        if name and sequence_parts:
            sequence = "".join(sequence_parts)
            outfile.write(f"{name}\t{sequence}\n")

parse_hg_tables()

import pandas as pd

def parse_bed12_line(line):
    """
    Parse a single BED12 line (string or list) into components.
    Removes specified prefixes from transcript names.
    Returns transcript name, block count, block sizes, and block starts.
    """
    if isinstance(line, str):
        fields = line.strip().split('\t')
    else:
        fields = line

    if len(fields) < 12:
        raise ValueError("Each BED12 line must have at least 12 columns")

    # Extract key fields
    name = fields[3]
    block_count = int(fields[9])
    block_sizes = list(map(int, fields[10].rstrip(',').split(',')))
    block_starts = list(map(int, fields[11].rstrip(',').split(',')))

    # Remove unwanted prefixes from transcript name
    prefixes_to_remove = [
        "mm39_ct_allisoforms_2790_",
        "CACNA1C_"
    ]
    for prefix in prefixes_to_remove:
        if name.startswith(prefix):
            name = name[len(prefix):]

    return name, block_count, block_sizes, block_starts

def get_splice_junctions(block_sizes):
    """
    Given a list of block sizes, return transcript-level indexes
    of block ends (splice junctions), excluding the last block end.
    """
    junctions = []
    current_pos = 0
    for size in block_sizes[:-1]:  # skip last block
        current_pos += size
        junctions.append(current_pos)
    return junctions


def compute_splice_junctions_from_bed12(bed12_input):
    """
    Compute splice junction positions (in transcript coordinates)
    for each transcript in a BED12 file or DataFrame.

    Parameters
    ----------
    bed12_input : str or pd.DataFrame
        Path to BED12 file, or a DataFrame with at least 12 BED columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ['transcript', 'splice_junctions']
    """
    if isinstance(bed12_input, str):
        df = pd.read_csv(bed12_input, sep='\t', comment='#', header=None)
    elif isinstance(bed12_input, pd.DataFrame):
        df = bed12_input.copy()
    else:
        raise TypeError("Input must be a BED12 file path or DataFrame")

    results = []
    for _, row in df.iterrows():
        name, block_count, block_sizes, block_starts = parse_bed12_line(row)
        junctions = get_splice_junctions(block_sizes)
        junctions_str = ','.join(map(str, junctions)) if junctions else 'None'
        results.append((name, junctions_str))

    result_df = pd.DataFrame(results, columns=['ID', 'splice_junctions'])
    return result_df

df = pd.read_csv('hgTables_parsed.txt', sep='\t')

results = []
for idx, row in df.iterrows():
    seq_id = row.iloc[0]     # first column
    cdna = row.iloc[1]       # second column

    cleaned = clean_seq(cdna)
    best_frame, longest_length, stop_index = find_longest_m_peptide(cleaned, translate_seq)
    best_frame = best_frame +1

    results.append({
        "ID": seq_id,
        "cDNA_length": len(cleaned),
        "best_frame": best_frame,
        "longest_length": longest_length,
        "stop_index": stop_index
    })

output_df = pd.DataFrame(results)

merged_df = output_df.merge(splicejunctions, on="ID", how="left")
def check_nmd(row):
    if pd.isna(row['splice_junctions']):
        return "Yes"  # No splice junctions to compare
    
    try:
        splice_sites = [int(x) for x in str(row['splice_junctions']).split(',') if x.strip().isdigit()]
    except ValueError:
        return "Yes"  # skip invalid entries
    
    stop_index = row['stop_index']
    
    # Check if any splice junction is within ±55 of stop_index
    for site in splice_sites:
        if abs(site - stop_index) <= 55:
            return "No"
    
    return "Yes"

# Apply to merged_df
merged_df["NMD"] = merged_df.apply(check_nmd, axis=1)

# Drop splice_junctions and print
final_df = merged_df.drop(columns=["splice_junctions"])
final_df.to_csv("output.tsv", sep="\t", index=False)

