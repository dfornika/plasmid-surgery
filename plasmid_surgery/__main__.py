import argparse
import csv
import json
import logging
import os

from textwrap import wrap

import plasmid_surgery.parsers as parsers
import plasmid_surgery.alignment as alignment

logger = logging.getLogger(__name__)

def setup_logging(log_level):
    """
    """
    numeric_log_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    logging.basicConfig(format='%(asctime)s :: %(levelname)s :: %(message)s', datefmt='%Y-%m-%dT%H:%M:%S', encoding='utf-8', level=numeric_log_level)

    return None


def parse_args():
    """
    """
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

    return args

def main():
    args = parse_args()
    setup_logging(args.log_level)
    output_dir = os.path.dirname(args.output)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    summary_output_dir = os.path.dirname(args.summary_output)
    if not os.path.exists(summary_output_dir):
        os.makedirs(summary_output_dir)

    summary_output = []

    plasmid_reconstruction_contigs = parsers.parse_fasta(args.plasmid_reconstruction)
    ref_plasmid = parsers.parse_fasta(args.ref_plasmid)[0]
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
        aligned_segments = alignment.align_contig_vs_ref(contig, ref_plasmid, args.tmpdir)
        
        if len(aligned_segments) == 0:
            logger.info(f"Zero alignments for contig: {contig_id} vs. ref_plasmid: {ref_plasmid_id}")
            continue
        elif len(aligned_segments) == 2 and aligned_segments[1]['is_supplementary']:
            logger.info(f"Supplementary alignment found for contig: {contig_id}")
            # TODO: Handle multiple alignments of contig vs ref
        
        for aligned_segment in aligned_segments:

            aligned_segment_qc_stats = alignment.alignment_qc(aligned_segment)
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
