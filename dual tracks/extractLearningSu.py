# -*- coding: utf-8 -*-
"""

To extract learning data in different data structures


"""

import os
import h5py
import re
import numpy as np
import pandas as pd
import su_region_align as align
import selectivity as zpy




trials = []
learningTrials = []
paths = []
ids = []
regions = []
unlabeledRecord = []

(regionL,non) = align.getRegionList(cmp=False)

for path in align.traverse(r"/home/zx/neupix/wyt/DataSum"):

    if not os.path.isfile(os.path.join(path, "events.hdf5")):
            print("missing one event file, path: ", path)
            continue

    print(path)

    # judge performance first
    os.chdir(path)
    trials=np.empty([0])
    with h5py.File(os.path.join(path, "events.hdf5"),'r') as fe:
            dset=fe['trials']
            trials=np.array(dset,dtype='int32')
    inWindow=align.judgePerformance(trials)
    if np.count_nonzero(inWindow == 2) >= 40:
        (bs_id, time_s, who) = align.get_bsid_duration_who(path)
        (mice_id, date, imec_no) = zpy.get_miceid_date_imecno(path)
        depthL = align.getTrackRegion(regionL, mice_id, date, imec_no, who)
        unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")
        wf = unitInfo[:].get('KSLabel')=='good'
        spkRate = unitInfo.iloc[:,7]
        try:
                spkRate = spkRate.str.replace(' spk/s','').astype(np.float64)
        except:
                spkRate = spkRate
        finally:
                goodSuIdx = unitInfo.iloc[:,[0]][wf & (spkRate > 1)].id
        # go detail into each unit
        for oneSuIdx in goodSuIdx:
                depth = unitInfo.loc[unitInfo['id'] == oneSuIdx, ['depth']].iat[0,0]          
                reg = align.matchDepth(depth, depthL, date, mice_id, imec_no, unlabeledRecord)
                learningTrials.append(inWindow)
                paths.append(path)
                ids.append(oneSuIdx)
                regions.append(reg)

os.chdir(r"/home/xd/data/learning")
f = h5py.File("learningSumFull.hdf5", "w")
string_dt = h5py.special_dtype(vlen=str)
d1 = f.create_dataset("/path",data=np.array(paths,dtype=object),dtype=string_dt)
d2 = f.create_dataset("/SUid", data=np.array(ids),dtype=int)
d3 = f.create_dataset("/reg", data=np.array(regions,dtype=object),dtype=string_dt)
f.close()


import scipy.io as scio
scio.savemat('learningTrialsFull.mat',{'learningTrials':learningTrials})
