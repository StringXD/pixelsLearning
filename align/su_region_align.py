# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:29:29 2019

@author: Libra

This script is responsible for generating the su_id2reg.csv which is further used by other statistics


"""

import os
import h5py
import sys
import re
import csv
import numpy as np
import pandas as pd
import fileutil as futil
import parsetree as parsetree

def imecNo2side(who_did, date, imecNo, mid): # assign probe # based on lab records
    if date == "191130" and mid == "49":
        return "L"

    if date == "191101" and mid == "26":
        return "R"

    if who_did == "HEM" and (int(date)) >= 191028:
        if imecNo == "1":
            return "R"
        elif imecNo == "0":
            return "L"
    else:
        if imecNo == "1":
            return "L"
        elif imecNo == "0":
            return "R"

    print("Error parsing imec No")
    return "X"


def combineSubRegion(r):
    if re.match("CA[13]", r):
        return r
    if re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r):
        g = re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r)
        return g.group(1)
    if r.startswith("COAp"):
        return "COAp"
    else:
        return r


def matchDepth(depth, depthL, date, mice_id, imec_no): # depth to brain-region-leaf
    if not depthL.empty:
        label = depthL.loc[
            (depthL["distance2tipLow"] - 100 <= depth) & (depth < depthL["distance2tipHigh"] + 140), # empirical 50% confidence interval
            ["acronym"],
        ]
        label = list(label['acronym'])
        labelparent = [combineSubRegion(subregion) for subregion in label]
        if len(set(labelparent)) == 1:
            return (label[0],None)
        else:
            ('unlabeled',None)

    unlabeledRecord=[date, mice_id, imec_no, depth]
    return ('unlabeled',unlabeledRecord)


def getTrackRegion(regionL, mice_id, date, imecNo): # all regions along a probe
    depthL = regionL.loc[
        (regionL["mouse_id"] == mice_id)
        & (regionL["implanting_date"] == date)
        & (regionL["imec"] == imecNo),
        ["acronym", "distance2tipLow", "distance2tipHigh"],
    ]
    return depthL


def getRegionList(): # file interface
    site_file = r"E:\prJ\neuropixels\learning\learningTracks.CSV"
    regionL = pd.read_csv(site_file).astype(
        {"mouse_id": "str", "implanting_date": "str", "imec": "str"}
    )[
        [
            "mouse_id",
            "implanting_date",
            "imec",
            "acronym",
            "distance2tipHigh",
            "distance2tipLow",
        ]
    ]

    return regionL


def align_onefolder(su_ids,path,regionL,offset): # processing one session
        (mice_id, date, imec_no) = futil.get_miceid_date_imecno(path) # meta data extracted from path and file tag
        depthL = getTrackRegion(regionL, mice_id, '20' + date, imec_no) # per probe brain regions
        su_region_corr = []
        unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), # for data cleaning
                               sep="\t",
                               usecols=['id','depth'],
                               index_col='id')
        unlabeled=[]

        ### sequentially match su ids with brain region tree
        for one_su in su_ids:
            depth = unitInfo.loc[one_su, ['depth']][0]
            (reg,one_unlabeled) = matchDepth(depth, depthL, date, mice_id, imec_no)
            if reg=='unlabeled':
                reg_depth=-1
                tree_str=['',]
            else:
                (regidx,reg_depth,tree_idces,tree_str)=parsetree.get_tree_path(reg)
            while len(tree_str)<6:
                tree_str.append('')
            su_region_corr.append([one_su+offset,reg_depth]+tree_str[:6])
            if one_unlabeled:
                unlabeled.append(one_unlabeled)
        return (su_region_corr,unlabeled)

### traverse all folder
def gen_align_files(cmp=False):
    regionL = getRegionList()
    unlabeledRecord = []

    for path in futil.traverse(r"E:\prJ\neuropixels\learning\dataSumLN"):
        su_region_corr = []
        print(path)
        with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr: # read back pre-cleaned data
            if not "SU_id" in ffr.keys():
                print("missing su_id key in path ", path)
                continue
            su_ids = np.array(ffr["SU_id"], dtype="double")[0]

        # Do not depend on imec number
        (su_region,unlabeled)=align_onefolder(su_ids,path,regionL,0)
        su_region_corr.append(su_region)
        unlabeledRecord.append(unlabeled)        

        su_region_corr=list(filter(None, su_region_corr))
        if su_region_corr: # result export to file
            tbl = pd.DataFrame(np.vstack(su_region_corr),
                               columns=['index', 'depth','d3','d4','d5','d6','d7','d8']).set_index('index')[:]
            tbl.to_csv(os.path.join(path, 'su_id2reg.csv'), header=True)
        # else:
            # breakpoint()

    unlabeledRecord=list(filter(None, unlabeledRecord))
    with open("unlabeled.csv", "w", newline="") as f: # Incomplete data, developers only.
        writer = csv.writer(f)
        writer.writerows(unlabeledRecord)