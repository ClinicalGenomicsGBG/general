#!/usr/bin/env python

"""Script for chunking fasta file sequences on their GC-content.

Takes one GC-content span, calculates sequences GC-content and prints into 'GC-content span'-file and 'remaining'-file.
This script is designed to be run as an external app in CLC.
"""

import re
import statistics
import tempfile
import sys
import subprocess
import shutil
from tabulate import tabulate
from Bio import SeqIO


class Context():
    def __init__(self, fasta_input, span_interest, span_output, remain_output, summary_output, length):
        self.fasta_input = fasta_input
        self.span_interest = span_interest
        self.span_output = span_output
        self.remain_output = remain_output
        self.length = int(length)
        self.summary_output = summary_output

        self.tmp_dir = tempfile.mkdtemp(prefix='gcsort')
        self.span_tmp = f'{self.tmp_dir}/span_tmp.fna'
        self.remain_tmp = f'{self.tmp_dir}/remain_tmp.fna'
        self.summary_tmp = f'{self.tmp_dir}/summary_tmp.rst'

        self.seq_records = self._read_fasta()  # Generator to save memory

        self.span_data = {self.span_interest: [],
                          'remaining': []}
        self.shorts = 0

    def _read_fasta(self):
        """Read fasta and create generator."""
        for record in SeqIO.parse(self.fasta_input, 'fasta'):
            yield record


def _verify_spans(span):
    """Verify the structure of the input argument for span of interest."""
    try:
        first, last = [float(i) for i in re.split('-', span)]
    except ValueError:
        sys.exit('Incorrect span: Only structure "10-20.5" allowed')

    if first == last:
        sys.exit('Incorrect span: Lower == Upper')

    if first > last:
        first, last = last, first

    # Should not be possible
    if first < 0 or last < 0:
        sys.exit('Incorrect percentage: < 0')

    if first > 100 or last > 100:
        sys.exit('Incorrect percentage: > 100')

    return (first/100, last/100)


def _verify_fasta_input(fasta):
    """Make sure supposed empty files are empty."""
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.seq:
            return fasta

    open(fasta, 'w').close()  # Blanks file
    return fasta


def parse_n_chunk(ctx):
    """Parse the fasta and write contigs to their respective outputs."""
    with open(ctx.span_tmp, 'w') as span_out, open(ctx.remain_tmp, 'w') as rem_out:
        for record in ctx.seq_records:
            gc = record.seq.upper().count('G') + record.seq.upper().count('C')
            gc_content = round(gc / len(record.seq), 4)

            if len(record.seq) < ctx.length:
                ctx.shorts += 1
                continue

            if ctx.span_interest[0] < gc_content <= ctx.span_interest[1]:
                SeqIO.write(record, span_out, 'fasta')
                ctx.span_data[ctx.span_interest].append(gc_content)
            else:
                SeqIO.write(record, rem_out, 'fasta')
                ctx.span_data['remaining'].append(gc_content)


def summarize(ctx):
    """Print short summary of GC-sort."""
    headers = ['Span', 'Number', 'Median GC']
    table_info = []

    with open(ctx.summary_tmp, 'w') as out:
        for span, gcs in ctx.span_data.items():
            if span != 'remaining':
                span = f'{span[0]} - {span[1]}'

            numb = len(gcs)
            if gcs:
                median_gc = statistics.median(gcs)
            else:
                median_gc = 0

            table_info.append([span, numb, median_gc])
        print(tabulate(table_info, headers=headers, tablefmt='rst'), file=out)
        print(f'Min length removed: {ctx.shorts}', file=out)


def reorganize(ctx):
    """Replace CLC expected outputs with result files."""
    subprocess.run(['mv', ctx.span_tmp, ctx.span_output])
    subprocess.run(['mv', ctx.remain_tmp, ctx.remain_output])
    subprocess.run(['mv', ctx.summary_tmp, ctx.summary_output])
    shutil.rmtree(ctx.tmp_dir)


def main():
    fasta_input, span_interest, span_output, remain_output, summary_output, length = sys.argv[1:]
    span_interest = _verify_spans(span_interest)
    fasta_input = _verify_fasta_input(fasta_input)
    ctx = Context(fasta_input, span_interest, span_output, remain_output, summary_output, length)

    parse_n_chunk(ctx)
    summarize(ctx)
    reorganize(ctx)


main()
