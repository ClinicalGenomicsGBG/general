#!/usr/bin/env python
"""Ion torrent server script for managing the medstore upload of Ion Torrent results."""

import os
import json
import sys
import subprocess


class Context():
    def __init__(self, result_path):
        self.result_path = result_path
        self.run_name = self.result_path.strip('/').split('/')[-1]
        self.basecall_path = f'{result_path}/basecaller_results'
        self.medstore_path = f'/results/medstore/{self.run_name}'


def upload(ctx):
    """Check if the run has already been converted and uploaded to medstore."""
    def parse_log(ctx):
        """Parse result log for relevant files."""
        basecaller_log = f'{ctx.basecall_path}/datasets_basecaller.json'

        with open(basecaller_log) as inp:
            jn = json.load(inp)

        good_bams = []
        for entry in jn['datasets']:
            if not entry['dataset_name'].startswith('none'):
                good_bams.append((entry['basecaller_bam'], entry['dataset_name']))

        return good_bams

    def on_medstore(ctx):
        """Check if already on medstore."""
        if os.path.isdir(ctx.medstore_path):
            return True
        return False

    def convert_and_move(ctx, good_bams):
        """Convert to fastq and move."""
        os.mkdir(ctx.medstore_path)

        for bam, name in good_bams:
            bam_path = f'{ctx.basecall_path}/{bam}'
            short_name = name.split('/')[0].replace(' ', '_')
            medstore_path = f'{ctx.medstore_path}/{bam}_{short_name}.fastq'
            subprocess.run(f'bedtools bamtofastq -i {bam_path} -fq {medstore_path}', shell=True)

    good_bams = parse_log(ctx)

    if not on_medstore(ctx):
        convert_and_move(ctx, good_bams)


def medstore_upload():
    result_path = sys.argv[1]
    ctx = Context(result_path)

    upload(ctx)


medstore_upload()
