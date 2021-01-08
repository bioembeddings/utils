#!/bin/python3

# The idea of this script is to take one input sequence set and randomly sample sequences from it,
# then run the pipeline on that random subset of sequences.

import argparse
import random
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from Bio import SeqIO
from bio_embeddings.utilities.config import read_config_file
from bio_embeddings.utilities.pipeline import execute_pipeline_from_config

logger = logging.getLogger(__name__)

# Args
parser = argparse.ArgumentParser(description='Embeds random subset of a sequence file.')

parser.add_argument('config_path', metavar='/path/to/pipeline_definition.yml', type=str, nargs=1,
                    help='The path to the config. For examples, see folder "parameter examples".')
arguments = parser.parse_args()
# Options
config_path = arguments.config_path[0]
config = read_config_file(config_path)
sequence_path = config['global']['sequences_file']
max_number_of_sequences = config['global'].get('max_number_of_sequences', 250)
max_len = config['global'].get('max_len', 100)
min_len = config['global'].get('min_len', 50)
#
filtered_sequences = list()
total_aa = 0
#

for seq_record in SeqIO.parse(sequence_path, "fasta"):
    if max_len > len(seq_record) > min_len:
        filtered_sequences.append(seq_record)

random_sample = random.sample(filtered_sequences, max_number_of_sequences)

for sequence in random_sample:
    total_aa += len(sequence)

logger.info(f"Total AA={total_aa}.")
logger.info(f"Total per-AA embedding size={4*total_aa*1024*3*pow(10, -6)}MB")

with TemporaryDirectory() as workdir:
    seq_path = Path(workdir).joinpath("sequences.fasta")

    SeqIO.write(random_sample, seq_path, "fasta")

    # Add sampled sequences
    config['global']['sequences_file'] = str(seq_path)

    logger.info("------ Starting pipeline execution...")
    execute_pipeline_from_config(config)
    logger.info("------ Finished pipeline execution.")