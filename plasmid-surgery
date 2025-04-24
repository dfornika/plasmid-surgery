#!/usr/bin/env python

import argparse
import csv
import json
import logging
import os
import subprocess
import re

from pathlib import Path
from textwrap import wrap

import pysam


def parse_fasta(fasta_path: Path):
    """
    Parse a fasta file into a list of dicts.

    :param fasta_path: Path to the fasta file
    :type fasta_path: Path
    :return: The parsed fasta file. List of dicts with keys: 'seq_id', 'description', 'seq', 'length'
    :rtype: list[dict]
    """
    seqs = []
    seq = {'seq': ""}
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if 'seq_id' in seq:
                    seq['length'] = len(seq['seq'])
                    seqs.append(seq)
                seq = {'seq': ""}
                defline = line.lstrip('>')
                seq_id = defline.split()[0]
                seq['seq_id'] = seq_id
                description = ' '.join(defline.lstrip(seq_id).split())
                seq['description'] = description
            else:
                seq['seq'] += line

    seq['length'] = len(seq['seq'])
    seqs.append(seq)

    return seqs


def revcomp(seq: str):
    """
    Reverse complement a DNA sequence (complements are: A -> T, C -> G, G -> C, T -> A)

    :param seq: The sequence to reverse complement.
    :type seq: str
    :return: The reverse complemented sequence.
    :rtype: str
    """
    reverse_complement_seq = ""
    lookup = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'a': 't',
        'c': 'g',
        'g': 'c',
        't': 'a',
        'N': 'N',
        'n': 'n'
    }
    reversed_seq = seq[::-1]
    for x in reversed_seq:
        c = lookup.get(x, None)
        if not c:
            c = 'N'
        reverse_complement_seq += c

    return reverse_complement_seq


def parse_cigar(cigar_str: str) -> list[tuple]:
    """
    Parse a CIGAR string into a list of (operation, length) tuples.
    Raises ValueError if the string is malformed.

    :param cigar_str:
    :type cigar_str:
    :return:
    :rtype:
    """
    cigar_regex = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = cigar_regex.findall(cigar_str)
    total_len = sum(len(m[0]) + 1 for m in matches)
    if total_len != len(cigar_str):
        raise ValueError(f"Invalid CIGAR string: {cigar_str}")

    cigar_operation_codes = {
        'M': 'alignment_match',
        'I': 'insertion_to_ref',
        'D': 'deletion_to_ref',
        'N': 'skip_ref',
        'S': 'soft_clip',
        'H': 'hard_clip',
        'P': 'padding',
        '=': 'sequence_match',
        'X': 'sequence_mismatch',
        'B': 'back'
    }
    
    parsed_cigar = [(cigar_operation_codes[op], int(length)) for length, op in matches]

    return parsed_cigar


def parse_sam(sam_path) -> list[dict]:
    """
    Parse a sam file into a list of dicts.
    """
    samfile = pysam.AlignmentFile(sam_path, "r")
    aligned_segments = []
    for r in samfile:
        contig_id = r.query_name
        ref_id = samfile.getrname(r.reference_id)  # name of reference
        ref_start = r.reference_start
        ref_end = r.reference_end
        query_start = r.query_alignment_start
        query_end = r.query_alignment_end
        query_alignment_length = r.query_alignment_length
        query_alignment_sequence = r.query_alignment_sequence
        is_reverse = r.is_reverse
        is_supplementary = r.is_supplementary
        is_mapped = r.is_mapped
        cigar_str = r.cigarstring
        cigar_parsed = []
        if cigar_str:
            cigar_parsed = parse_cigar(cigar_str)
        tags = dict(r.tags)
        cg = tags.get("cg", None)
        edit_distance = tags.get("NM", None)

        matches = 0
        if r.cigartuples:
            matches = sum(length for (op, length) in r.cigartuples if op == 7)  # '=' in CIGAR
        mismatches = 0
        if r.cigartuples:
            mismatches = sum(length for (op, length) in r.cigartuples if op == 8)  # 'X' in CIGAR
        percent_identity = 0.0
        if (matches + mismatches) > 0:
            percent_identity = round(100 * matches / (matches + mismatches), 3)
        
        aligned_segment = {
            'contig_id': contig_id,
            'ref_id': ref_id,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'query_start': query_start,
            'query_end': query_end,
            'query_alignment_length': query_alignment_length,
            'query_alignment_sequence': query_alignment_sequence,
            'is_reverse': is_reverse,
            'is_supplementary': is_supplementary,
            'is_mapped': is_mapped,
            'cigar_str': cigar_str,
            'cigar_parsed': cigar_parsed,
            'tags': tags,
            'edit_distance': edit_distance,
            'percent_identity': percent_identity,
        }
        aligned_segments.append(aligned_segment)

    # print(json.dumps(aligned_segments, indent=2))
        
    return aligned_segments


def minimap2_align_contig_vs_ref(contig: dict, ref: dict, tmpdir: Path) -> list[dict]:
    """
    Align contig vs ref using minimap2

    :param contig:
    :type contig: dict
    :param ref:
    :type ref: dict
    :param tmpdir:
    :type tmpdir: Path
    :return:
    :rtype: list[dict]
    """
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    contig_id = contig['seq_id']
    ref_id = ref['seq_id']
    contig_path = os.path.join(tmpdir, f"{contig_id}.fa")
    
    with open(contig_path, 'w') as f:
        f.write(f">{contig['seq_id']}\n")
        f.write(f"{contig['seq']}\n")

    ref_path = os.path.join(tmpdir, f"{ref_id}.fa")

    if not os.path.exists(ref_path):
        with open(ref_path, 'w') as f:
            f.write(f">{ref['seq_id']} [circular=true]\n")
            f.write(f"{ref['seq']}\n")
    alignment_path = os.path.join(tmpdir, f"{contig_id}_vs_{ref_id}.sam")
    
    minimap2_cmd = [
        'minimap2',
        '-a',
        '-cx',
        'asm10',
        '--secondary=no',
        '--eqx',
        '-o', alignment_path,
        ref_path,
        contig_path,
    ]

    try:
        subprocess.run(minimap2_cmd, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e)
        exit(-1)

    aligned_segments = parse_sam(alignment_path)
    for aligned_segment in aligned_segments:
        aligned_segment['query_full_length'] = len(contig['seq'])
        aligned_segment['ref_full_length'] = len(ref['seq'])

    return aligned_segments


def alignment_qc(aln, min_identity=80.0, min_aln_length=50, max_gap_frac=0.8):
    """
    Determine if alignment passes filters.
    
    """
    gap_bases = sum(length for op, length in aln['cigar_parsed'] if op in ('insertion_to_ref', 'deletion_to_ref'))
    gap_frac = gap_bases / aln['query_alignment_length']

    alignment_qc = {
        'qc_stats': {
            'min_identity': {
                'min_value': min_identity,
                'value': aln['percent_identity'],
            },
            'aln_length': {
                'min_value': min_aln_length,
                'value': aln['query_alignment_length'],
            },
            'gap_frac': {
                'max_value': max_gap_frac,
                'value': gap_frac,
            },
        },
        'overall_qc_pass': True,
    }
    
    for qc_stat, qc_values in alignment_qc['qc_stats'].items():
        actual_value = qc_values['value']
        if 'min_value' in qc_values:
            if actual_value < qc_values['min_value']:
                alignment_qc['overall_qc_pass'] = False
        elif 'max_value' in qc_values:
            if actual_value > qc_values['max_value']:
                alignment_qc['overall_qc_pass'] = False
    
    return alignment_qc


def main(args):
    logger = logging.getLogger(__name__)
    numeric_log_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError(f"Invalid log level: {args.log_level}")

    logging.basicConfig(format='%(asctime)s :: %(levelname)s :: %(message)s', datefmt='%Y-%m-%dT%H:%M:%S', encoding='utf-8', level=numeric_log_level)
    
    output_dir = os.path.dirname(args.output)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    summary_output_dir = os.path.dirname(args.summary_output)
    if not os.path.exists(summary_output_dir):
        os.makedirs(summary_output_dir)
        
    summary_output = []
    
    plasmid_reconstruction_contigs = parse_fasta(args.plasmid_reconstruction)
    ref_plasmid = parse_fasta(args.ref_plasmid)[0]
    ref_plasmid_length = len(ref_plasmid['seq'])
    ref_plasmid_id = ref_plasmid['seq_id']

    output_seq = 'N' * ref_plasmid_length
    output_seq_id = None
    if args.output_seq_id:
        output_seq_id = args.output_seq_id
    else:
        plasmid_reconstruction_basename = os.path.basename(args.plasmid_reconstruction)
        output_seq_id = plasmid_reconstruction_basename.split('.')[0]


    aligned_contigs = []
    failed_alignments = []
    for contig in plasmid_reconstruction_contigs:
        contig_id = contig['seq_id']
        contig['reverse_complement'] = False
        logger.info(f"Aligning contig: {contig_id} vs. ref plasmid: {ref_plasmid_id}")
        aligned_segments = minimap2_align_contig_vs_ref(contig, ref_plasmid, args.tmpdir)
        
        if len(aligned_segments) == 0:
            logger.info(f"Zero alignments for contig: {contig_id} vs. ref_plasmid: {ref_plasmid_id}")
            continue
        elif len(aligned_segments) == 2 and aligned_segments[1]['is_supplementary']:
            logger.info(f"Supplementary alignment found for contig: {contig_id}")
            # TODO: Handle multiple alignments of contig vs ref
        
        for aligned_segment in aligned_segments:

            aligned_segment_qc_stats = alignment_qc(aligned_segment)
            passes_alignment_filters = aligned_segment_qc_stats['overall_qc_pass']
            if not(passes_alignment_filters):
                logger.info(f"Alignment failed quality filters: {contig_id}")
                logger.info("Alignment QC stats: " + json.dumps(aligned_segment_qc_stats['qc_stats']))
                aligned_segment['included_in_output'] = False
                failed_alignments.append(aligned_segment)
                continue

            aligned_contigs.append(aligned_segment)


    max_replacement_idx = ref_plasmid_length + 1
    for aligned_contig in aligned_contigs:
        aligned_contig['included_in_output'] = False
        replacement_start = aligned_contig['ref_start']
        replacement_end = aligned_contig['ref_end']
        if replacement_start <= max_replacement_idx and replacement_end <= max_replacement_idx:
            aligned_contig['included_in_output'] = True
            output_seq = output_seq[:replacement_start] + aligned_contig['query_alignment_sequence'] + output_seq[replacement_end:]

    # check outputs
    #for aligned_contig in aligned_contigs:        
    #    aligned_contig.pop('ref')
    #    aligned_contig.pop('query')
    #    print(json.dumps(aligned_contig, indent=2))

    for aligned_contig in aligned_contigs + failed_alignments:
        contig_id = aligned_contig.get('contig_id', None)
        reverse_complement = aligned_contig.get('is_reverse', None)
        is_mapped = aligned_contig.get('is_mapped', None)
        ref_id = ref_plasmid_id
        alignment_start = aligned_contig.get('ref_start', None)
        alignment_end = aligned_contig.get('ref_end', None)
        alignment_length = aligned_contig.get('query_alignment_length', None)
        if not is_mapped:
            alignment_length = 'NA'
        percent_identity = aligned_contig.get('percent_identity', None)
        included_in_output = aligned_contig.get('included_in_output', None)

        FAKE_ALIGNMENT_START = 1_000_000_000_000
        if not is_mapped:
            alignment_start = FAKE_ALIGNMENT_START
            alignment_end = 'NA'
            alignment_length = 'NA'
            percent_identity = 'NA'

        summary_output_row = {
            'contig_id': contig_id,
            'ref_id': ref_id,
            'ref_length': ref_plasmid_length,
            'reverse_complement': reverse_complement,
            'is_mapped': is_mapped,
            'alignment_start': alignment_start,
            'alignment_end': alignment_end,
            'alignment_length': alignment_length,
            'percent_identity': percent_identity,
            'included_in_output': included_in_output,
        }
        summary_output.append(summary_output_row)

    summary_output_fieldnames = [
        'contig_id',
        'ref_id',
        'ref_length',
        'reverse_complement',
        'is_mapped',
        'alignment_start',
        'alignment_end',
        'alignment_length',
        'percent_identity',
        'included_in_output',
    ]
    
    summary_output_sorted = sorted(summary_output, key=lambda x: x.get('alignment_start', 0))
    with open(args.summary_output, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=summary_output_fieldnames, delimiter='\t', extrasaction='ignore', quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for row in summary_output_sorted:
            if row['alignment_start'] == FAKE_ALIGNMENT_START:
                row['alignment_start'] = 'NA'
            writer.writerow(row)

    output_seq_wrapped = wrap(output_seq, args.output_bases_per_line)
    with open(args.output, 'w') as f:
        f.write(f">{output_seq_id}\n")
        for line in output_seq_wrapped:
            f.write(f"{line}\n")
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--plasmid-reconstruction', help="Plasmid reconstruction to be repaired.")
    parser.add_argument('--ref-plasmid', help="Reference plasmid to be used for plasmid surgery.")
    parser.add_argument('--tmpdir', default="./plasmid-surgery_tmp", help="temp dir to use for minimap2 alignments (default: './plasmid-surgery_tmp')")
    parser.add_argument('--summary-output', default='./plasmid-surgery-summary.tsv', help="Summary file for plasmid surgery process. (default: './plasmid-surgery-summary.tsv')")
    parser.add_argument('--output', default='./plasmid-surgery-output.fa', help="Path to output file for repaired plasmid reconstruction (default: './plasmid-surgery-output.fa')")
    parser.add_argument('--output-seq-id', help="Sequence ID to use in the output fasta header. If not supplied, uses plasmid reconstruction filename (minus file extension)")
    parser.add_argument('--output-bases-per-line', default=70, help="Number of bases per line in the output fasta (default: 70)")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    args = parser.parse_args()
    main(args)
