#!/usr/bin/env python
"""
Script for parsing the best hit references from a mapping stage of several genes with several variants.
"""

import csv
import sys
import itertools

from Bio import SeqIO


class MissingFastaError(object):
    pass


def parse_mapping(mapping_results):
    """Parse a vcf for the best references according to total read count."""
    with open(mapping_results) as inp:
        csv_info = csv.DictReader(inp, delimiter=',')
        rows = [row for row in csv_info]

    genegroups = [list(g) for _, g in itertools.groupby(rows, lambda x: x['Name'].split('_')[0])]
    best_mapped_genes = [max(group, key=lambda x: int(x['Total read count'])) for group in genegroups]
    return best_mapped_genes


def fetch_sequence_records(best_mapped_genes, fasta_path):
    """Fetch the SeqIO record for the best reference from references fasta file."""
    headers = [gene['Name'].replace(' mapping', '') for gene in best_mapped_genes]

    records = []
    for record in SeqIO.parse(fasta_path, 'fasta'):
        if record.id in headers:
            records.append(record)

    if len(records) != len(best_mapped_genes):
        raise MissingFastaError('One or more sequence headers were not found in the given fasta file.')
    return records


def write_output(gene_sequence_records, output):
    """Write the SeqIO record for the best reference to file."""
    with open(output, 'w') as out:
        for record in gene_sequence_records:
            SeqIO.write(record, out, 'fasta')


def main():
    mapping_results, fasta_path, fasta_output = sys.argv[1:]

    best_mapped_genes = parse_mapping(mapping_results)
    gene_sequence_records = fetch_sequence_records(best_mapped_genes, fasta_path)
    write_output(gene_sequence_records, fasta_output)


if __name__ == '__main__':
    main()
