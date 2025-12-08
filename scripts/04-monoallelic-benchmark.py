import os
SCRIPT_PATH   = os.path.dirname(__file__)
PROJROOT      = os.path.join(SCRIPT_PATH, '..')
TUNE_DATA_DIR = os.path.join(PROJROOT, 'data', 'model-tuning')

TUNE_RES_DIR = os.path.join(PROJROOT, 'results', '01-ablation-study')

import sys
sys.path.append(PROJROOT)

import argparse
parser = argparse.ArgumentParser()

import numpy  as np
import pandas as pd

from promiselib.model import PROMISE
from promiselib.Utils import readPickle
from promiselib.RunUtils import readPickle, model_load

MAX_PROTEIN_LEN = 3000
fbase = f"p{MAX_PROTEIN_LEN}"
pkl_file  = os.path.join(TUNE_DATA_DIR, f"{fbase}.pkl")
data = readPickle(pkl_file)

tbl_pkl_file  = os.path.join(TUNE_DATA_DIR, f"{fbase}-tables.pkl")
dat_tbl = readPickle(tbl_pkl_file)


BASE_TDIR = os.path.join(TUNE_RES_DIR, fbase)


TDIRS = [os.path.join(BASE_TDIR, d, 'seed-00') for d in os.listdir(BASE_TDIR)]

sp_ord = ['train', 'val', 'test']
for TDIR in TDIRS:
    fs = {f: os.path.join(TDIR, f) for f in os.listdir(TDIR)}
    mod_f = fs['model.keras']
    if os.path.exists(mod_f):
        mod = model_load(PROMISE, mod_f)
        res_lst = []
        for sp_name in sp_ord:
            inp, out = data[sp_name]
            out      = out[0][:,0]
            tbl      = dat_tbl[sp_name]
            if not all(tbl.loc[:,"label"] == out):
                raise ValueError("Mismatch between inputs")

            pred = mod.predict(inp)
            res_df = tbl.copy()
            res_df.loc[:,"pred"] = pred[:,0]
            res_lst.append(res_df)

        res_df = pd.concat(res_lst, axis = 0)
        res_f = os.path.join(TDIR, 'all-pred.tsv')
        res_df.to_csv(res_f, sep = '\t', index = False)
    else:
        print(f"Unable to find model file at: [mod_f]")
