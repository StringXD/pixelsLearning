# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra


@author: Libra


"""



import sys
import pickle

sys.path.insert(0, r"E:\prJ\neuropixels\learning\pixels-masterUpdate")

 
from per_sec_stats.prepare_data import prepare_data
from per_sec_stats.export import exporth5py
sys.path.append(r"E:\prJ\neuropixels\learning\pixels-masterUpdate\align")
from align import su_region_align


def gen_align_files():
    '''
    Generate su_id2reg csv file
    '''
    su_region_align.gen_align_files()

def gen_selectivity_stats(delay, debug = False, denovo = True):
    '''
    Generate per SU selectivity and brain region tree file
    '''
    if denovo:
        (dict_stats,error_files)=prepare_data(delay = delay, debug = debug)
        pickle.dump(dict_stats,open(f'per_sec_sel_{delay}.p','wb'))
        #TODO ^^^^ deprecated, for debug only
    else:
        dict_stats = pickle.load(open(f'per_sec_sel_{delay}.p','rb'))

    exporth5py(dict_stats) # transient_{delay}.hdf5
    return error_files


if __name__ == "__main__":
    delay = 6
    error_files=gen_selectivity_stats(delay, debug=False)