import tensorflow as tf
import tensorflow_hub as hub
import joblib
import gzip
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import time

import logging
import os
import sys

from Modules.Enformer import *
from Modules.FastaExt import *
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--fasta_file",
                    help="Relative path of the sequence fasta file (.fasta)", type=str, required=True)
parser.add_argument("--cell_type",
                    help="Cell type for experiment", type=str, required=True)
parser.add_argument("--experiment_name",
                    help="Name for experiment", type=str, required=True)
parser.add_argument("--indicies",
                    help="indicies for that cell type", type=int, nargs = '+', required=True)


hparams, _ = parser.parse_known_args()

start_time = time.time()
path = "/project/st-cdeboer-1/iluthra/enformer/test_enformer/trained_model/"
#model = EnformerScoreVariantsRaw(path) #run the model and store it in the variable
#
model = Enformer(path)

SEQUENCE_LENGTH = 393216

## read in fasta file

prefix = str("/scratch/st-cdeboer-1/iluthra/randomDNA/" + hparams.cell_type + '/')
input_fasta = open('/project/st-cdeboer-1/iluthra/enformer_random_DNA/enformer_data/' + hparams.fasta_file, 'r')

data = []

for line in input_fasta:
    if line.startswith('>'):
        continue
    else:
        data.append(line)


all_predictions = np.empty([896, len(hparams.indicies), 0])

for i in range(0, len(data)):
  #print(data[i])
  sequence_one_hot = one_hot_encode(data[i].strip())
  predictions = model.predict_on_batch(sequence_one_hot[np.newaxis])['human'][0]
  print(predictions.shape)
  predictions_current = predictions[:, hparams.indicies]
  print(predictions_current.shape)
  predictions_current = predictions_current[:,:, np.newaxis]
  print(predictions_current.shape)
  all_predictions = np.append(all_predictions, predictions_current, axis = 2)


  print("Predictions Done %d!" % (i))

##add if size != 896xn call error
print(time.time() - start_time)
print(all_predictions.shape)
np.save(prefix+hparams.experiment_name, all_predictions)
