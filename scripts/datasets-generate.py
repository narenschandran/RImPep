import os
import sys
import pickle
import numpy as np
import pandas as pd


projroot = os.path.join(os.path.dirname(__file__), '..')
prereq_dir = os.path.join(projroot, 'prereq')

sys.path.append(projroot)

from promiselib.Utils import readPickle, savePickle, pad_arr

model_inp_dir = os.path.join(projroot, 'data', 'model-tuning')

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

inp_dir = os.path.join(projroot, 'data', 'model-inputs')
sp_ds = {d: os.path.join(inp_dir, d) for d in os.listdir(inp_dir)}

seq_f = os.path.join(seq_cache_dir, 'human-proteome-canonical.seq.pkl')

# When I first started developing the code to collate
# the different datapoints for each split, I used
# os.listdir to fetch all the folders within each split.
# Turns out, the order in which os.listdir fetches the
# folders seems to be dependent on the specific machine/os
# that I run it on. In order to ensure that the tuning data
# is exactly identical across everyone's runs, I've fixed
# the sample order based on the order in the first machine
# where I developed this code.
sample_tune_f = os.path.join(prereq_dir, 'sample-tuning-shuffle.txt')
sample_tune_ord = {}
with open(sample_tune_f, 'r') as samp_fcon:
    for raw_samp_line in samp_fcon:
        samp_line = raw_samp_line.strip()
        sp_name, samps = samp_line.split(":")
        sample_tune_ord[sp_name] = samps.split(",")


seq_dct = readPickle(seq_f)

split_lst = {split_name: [] for split_name in sample_tune_ord.keys()}

out_dat = {}
for split_name in sample_tune_ord.keys():
    base_d = sp_ds[split_name]
    ds = {d: os.path.join(base_d, d) for d in os.listdir(base_d)}

    tpm_lst = []
    seq_lst = []
    hla_lst = []
    lab_lst = []
    for sample_id in sample_tune_ord[split_name]:
        d = ds[sample_id]
        print(sample_id)

        tpm_f = os.path.join(d, 'tpm.tsv')
        dat = pd.read_csv(tpm_f, sep = '\t')


        al_f  = os.path.join(d, 'alleles.txt')
        with open(al_f, 'r') as fcon:
            alleles = [line.strip() for line in fcon]


        for allele in alleles:
            hla_dct = readPickle(hla_pkls[allele])
            tpm_lst.append(np.reshape(dat.tpm.values, (dat.shape[0], 1)))
            seq_lst.append(np.asarray([seq_dct[gene] for gene in dat.gene]))
            hla = []
            for gene in dat.gene:
                hla1 = [float(x) for x in pad_arr(list(hla_dct[gene][:,0]), 3000)]
                hla2 = [float(x) for x in pad_arr(list(hla_dct[gene][:,1]), 3000)]
                hla.append(np.asarray([hla1, hla2]))
            hla_lst.append(np.transpose(np.asarray(hla), (0, 2, 1)))
            lab_lst.append(np.reshape(dat.label.values, (dat.shape[0], 1)))
            tmp_datf = dat.copy()
            tmp_datf.loc[:,"allele"] = allele
            tmp_datf.loc[:,"split_name"] = split_name
            split_lst[split_name].append(tmp_datf.copy())
            del tmp_datf

    out_dat[split_name] = [
        [
            np.concatenate(tpm_lst, axis = 0),
            np.concatenate(seq_lst, axis = 0),
            np.concatenate(hla_lst, axis = 0)
        ],
        [np.concatenate(lab_lst, axis = 0)]
    ]
    del tpm_lst
    del seq_lst
    del hla_lst
    del lab_lst



outf = os.path.join(model_inp_dir, 'p3000.pkl')
savePickle(out_dat, outf)


split_dat = {split_name: pd.concat(sp_dat, axis = 0) for split_name, sp_dat in split_lst.items()}
tbl_outf = os.path.join(model_inp_dir, 'p3000-tables.pkl')
savePickle(split_dat, tbl_outf)
