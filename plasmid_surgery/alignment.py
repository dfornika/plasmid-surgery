import os
import logging
import subprocess

from pathlib import Path

import plasmid_surgery.parsers as parsers

logger = logging.getLogger(__name__)

def align_contig_vs_ref(contig: dict, ref: dict, tmpdir: Path) -> list[dict]:
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

    aligned_segments = parsers.parse_sam(alignment_path)
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
