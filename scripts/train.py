import os
SCRIPT_PATH   = os.path.dirname(__file__)
PROJROOT      = os.path.join(SCRIPT_PATH, '..')
TUNE_DATA_DIR = os.path.join(PROJROOT, 'data', 'model-tuning')

import sys
sys.path.append(PROJROOT)

import argparse
parser = argparse.ArgumentParser()

#------------------------------------------------#
#         General configuration options          #
#------------------------------------------------#

parser.add_argument(
    "-s", "--seed",
    help = "Seed value for training model",
    required = True, type = int
)

parser.add_argument(
    "-f", "--force",
    help = "Force model training even if prior runs exist",
    action = "store_true",
)

parser.add_argument(
    "-d", "--output_dir",
    help = "Base output directory",
    required = False, type = str, default = 'tmp/model-tuning'
)

parser.add_argument(
    "-b", "--batch_size",
    help = "Batch size for model tuning",
    required = False, type = int, default = 256
)

#------------------------------------------------#
#          Model configuration options           #
#------------------------------------------------#
parser.add_argument(
    "-l", "--max_prot_len",
    help = "Maximum protein length",
    required = True, type = int
)

parser.add_argument(
    "-w", "--window_size",
    help = "Window size for hierarchical LSTM",
    required = True, type = int
)

# Module activation and associated hyperparameters
parser.add_argument(
    "-T", "--use_tpm",
    help = "Train model with TPM module",
    action = "store_true",
)

parser.add_argument(
    "-C", "--seq_annconf",
    required = False, default = None
)
parser.add_argument(
    "-S", "--seq_embdim",
    required = False, type = int, default = 0
)


parser.add_argument(
    "-D", "--seq_lstmdim",
    required = False, type = int, default = 0
)


parser.add_argument(
    "-E", "--epi_lstmdim",
    help = "Train model with epitope scanning module (specify hLSTM embdim here to activate)",
    required = False, type = int, default = 0
)

parser.add_argument(
    "-F", "--epi_annconf",
    required = False, default = None
)

parser.add_argument(
    "-O", "--prob_annconf",
    help = "Probability computation output layer",
    required = False, type = str,
    default = 'gelu64_sigmoid1',
)

args = parser.parse_args()


import numpy  as np
import pandas as pd

from promiselib.model import PROMISE
from promiselib.TrainUtils import model_train
from promiselib.Utils import readPickle

pkl_fname = 'p3000.pkl'
pkl_file  = os.path.join(TUNE_DATA_DIR, pkl_fname)
data = readPickle(pkl_file)

params = {
    'max_prot_size'    : args.max_prot_len,
    'wsize'            : args.window_size,
    'use_tpm'          : args.use_tpm,
    'seq_embdim'       : args.seq_embdim,
    'seq_lstmdim'      : args.seq_lstmdim,
    'seq_annconf'      : args.seq_annconf,
    'epi_lstmdim'      : args.epi_lstmdim,
    'epi_annconf'      : args.epi_annconf,
    'prob_annconf'     : args.prob_annconf,
}

print("#----------------#")
print("|Model parameters|")
print("#----------------#")
for k, v in params.items():
    print(f"{k}: {v}")
print("#----------------#")

hist, perf, mod = model_train(
    modfn          = PROMISE,
    params         = params,
    data           = data,
    seed           = args.seed,
    base_outdir    = args.output_dir,
    epochs         = 500,
    verbose        = 1,
    early_stopping = True,
    checkpoints    = True,
    tensorboard    = True,
    force          = args.force,
    batch_size     = args.batch_size
)
