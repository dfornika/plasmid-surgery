# Plasmid Surgery

Tools like [mob-suite](https://github.com/phac-nml/mob-suite) produce plasmid reconstructions that consist of a set of contigs
without information about the order or orientation of the contigs in the complete plasmid.

This tool takes such a plasmid reconstruction along with a suitable reference plasmid, and attempts to identify the order,
orientation and (approximate) location of each contig in the plasmid reconstruction, based on where it best aligns on the
reference plasmid.

## Usage

```
usage: plasmid-surgery.py [-h] [--plasmid-reconstruction PLASMID_RECONSTRUCTION] [--ref-plasmid REF_PLASMID] [--tmpdir TMPDIR] [--summary-output SUMMARY_OUTPUT] [--output OUTPUT] [--output-seq-id OUTPUT_SEQ_ID]
                          [--output-bases-per-line OUTPUT_BASES_PER_LINE] [--log-level LOG_LEVEL]

options:
  -h, --help            show this help message and exit
  --plasmid-reconstruction PLASMID_RECONSTRUCTION
                        Plasmid reconstruction to be repaired.
  --ref-plasmid REF_PLASMID
                        Reference plasmid to be used for plasmid surgery.
  --tmpdir TMPDIR       temp dir to use for minimap2 alignments (default: './plasmid-surgery_tmp')
  --summary-output SUMMARY_OUTPUT
                        Summary file for plasmid surgery process. (default: './plasmid-surgery-summary.tsv')
  --output OUTPUT       Path to output file for repaired plasmid reconstruction (default: './plasmid-surgery-output.fa')
  --output-seq-id OUTPUT_SEQ_ID
                        Sequence ID to use in the output fasta header. If not supplied, uses plasmid reconstruction filename (minus file extension)
  --output-bases-per-line OUTPUT_BASES_PER_LINE
                        Number of bases per line in the output fasta (default: 70)
  --log-level LOG_LEVEL
                        Log level (default: 'info')
```

Example:

```
plasmid-surgery \
  --plasmid-reconstruction your_plasmid_reconstruction.fa \
  --ref-plasmid your_ref_plasmid.fa \
  --output fixed_plasmid_reconstruction.fa \
  --summary-output plasmid-surgery-summary.tsv 
```

## Outputs

The plasmid surgery summary tsv file includes these fields:

| field                | description                                                             |
|----------------------|-------------------------------------------------------------------------|
| `contig_id`          | Contig ID from plasmid reconstruction                                   |
| `ref_id`             | Ref plasmid ID                                                          |
| `ref_length`         | Length of ref plasmid                                                   |
| `reverse_complement` | Aligned sequence is reverse-complement of original contig seq           |
| `is_mapped`          | Contig is mapped somewhere against the ref plasmid                      |
| `alignment_start`    | Location on the ref plasmid where the contig alignment starts           |
| `alignment_end`      | Location on the ref plasmid where the contig alignment ends             |
| `alignment_length`   | Length of the alignment                                                 |
| `percent_identity`   | Percent identity of the alignment of the contig against the ref plasmid |
| `included_in_output` | The contig is included in the output complete plasmid seq               |

