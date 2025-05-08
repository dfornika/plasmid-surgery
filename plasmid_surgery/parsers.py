import pysam
import logging
import re

from pathlib import Path

logger = logging.getLogger(__name__)

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
