import os
SCRIPT_PATH   = os.path.dirname(__file__)
PROJROOT      = os.path.join(SCRIPT_PATH, '..')
TUNE_DATA_DIR = os.path.join(PROJROOT, 'data', 'model-tuning')
TUNE_RES_DIR  = os.path.join(PROJROOT, 'results', '01-ablation-study')
DISSECT_DIR   = os.path.join(PROJROOT, 'results', '02-model-dissection')

import sys
sys.path.append(PROJROOT)


import sys
import pickle
import numpy as np
import pandas as pd


from promiselib.Utils import readPickle, savePickle, pad_arr
from promiselib.RunUtils import model_load
from promiselib.model import PROMISE

#---------------------------------------#
# Data setup

model_inp_dir = os.path.join(PROJROOT, 'data', 'model-tuning')

cache_dir = os.path.join(model_inp_dir, 'cache')
seq_cache_dir = os.path.join(cache_dir, 'seq')
scan_cache_dir = os.path.join(cache_dir, 'scan')

def allele_clean(y):
    x = os.path.splitext(os.path.basename(y))[0]
    y = os.path.splitext(x)[0].split('-')
    v = y[0] + '-' + y[1] + '*' + y[2] + ':' + y[3]
    return v

hla_pkl_dir = scan_cache_dir
hla_pkls   = {allele_clean(f) : os.path.join(hla_pkl_dir, f) for f in os.listdir(hla_pkl_dir) if f[-4:] == ".pkl"}

seq_f = os.path.join(seq_cache_dir, 'human-proteome-canonical.seq.pkl')
seq_dct = readPickle(seq_f)

hla_dct = readPickle(hla_pkls['HLA-A*02:01'])

genes = sorted(list(set(hla_dct.keys()).intersection(set(seq_dct.keys()))))
dat = pd.DataFrame({'gene' : genes})

def inps_from_dat(datf):

    inps = [
        np.reshape(datf.tpm.values, (datf.shape[0], 1)),
        np.concatenate(
            [np.asarray([seq_dct[gene] for gene in datf.gene])],
        axis = 0)]

    hla = []
    for gene in datf.gene:
        hla1 = [float(x) for x in pad_arr(list(hla_dct[gene][:,0]), 3000)]
        hla2 = [float(x) for x in pad_arr(list(hla_dct[gene][:,1]), 3000)]
        hla.append(np.asarray([hla1, hla2]))

    inps.append(np.transpose(np.asarray(hla), (0, 2, 1)))
    return inps

base_odir = os.path.join(PROJROOT, 'results', '02-model-dissection')
if not os.path.exists(base_odir):
    os.makedirs(base_odir)


odir1 = os.path.join(base_odir, 'tpm')
if not os.path.exists(odir1):
    os.makedirs(odir1)

tpm_datf = dat.copy()
tpm_datf.loc[:,"tpm"] = np.linspace(0.0, 1000.0, tpm_datf.shape[0])
savePickle(tpm_datf, os.path.join(odir1, 'input-table.pkl'))

tpm_inps = inps_from_dat(tpm_datf)
savePickle(tpm_inps, os.path.join(odir1, 'inputs.pkl'))


odir2 = os.path.join(base_odir, 'tpm-seq')
if not os.path.exists(odir2):
    os.makedirs(odir2)

tseq_datf = dat.copy()
tseq_datf.loc[:,"tpm"] = 100.0
savePickle(tseq_datf, os.path.join(odir2, 'input-table.pkl'))

tseq_inps = inps_from_dat(tseq_datf)
savePickle(tseq_inps, os.path.join(odir2, 'inputs.pkl'))

#---------------------------------------#
# Prediction 

MAX_PROTEIN_LEN = 3000
fbase = f"p{MAX_PROTEIN_LEN}"
BASE_TDIR = os.path.join(TUNE_RES_DIR, fbase)
TDIRS0    = [os.path.join(BASE_TDIR, d, 'seed-00') for d in os.listdir(BASE_TDIR)]


# We'll only be using the TPM only and
# Seq only models for dissection
TDIRS = {
    [TDIR for TDIR in TDIRS0 if 'txx' in TDIR][0] : odir1,
    [TDIR for TDIR in TDIRS0 if 'tsx' in TDIR][0] : odir2
}

DATFS = {
    [TDIR for TDIR in TDIRS0 if 'txx' in TDIR][0] : tpm_datf,
    [TDIR for TDIR in TDIRS0 if 'tsx' in TDIR][0] : tseq_datf,
}

for TDIR, pred_odir in TDIRS.items():
    mname = os.path.basename(os.path.dirname(TDIR))
    fs = {f: os.path.join(TDIR, f) for f in os.listdir(TDIR)}
    mod_f = fs['model.keras']
    if os.path.exists(mod_f):
        mod = model_load(PROMISE, mod_f)
        res_lst = []
        res_df   = DATFS[TDIR]
        pred = mod.predict(readPickle(os.path.join(pred_odir, 'inputs.pkl')))
        res_df.loc[:,"pred"] = pred[:,0]
        outf = os.path.join(pred_odir, "pred.tsv")
        res_df.to_csv(outf, sep = '\t', index = False)
    else:
        print(f"Unable to find model file at: [mod_f]")
