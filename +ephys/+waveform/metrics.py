# Copyright © 2019. Allen Institute.  All rights reserved.
# params ref: https://github.com/Julie-Fabre/bombcell/wiki/Detailed-overview-of-quality-metrics#number-of-waveform-peaks


import os
import numpy as np
import pandas as pd
from collections import OrderedDict
import math
import warnings

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score

import scipy.io as scio
from scipy.spatial.distance import cdist
from scipy.stats import chi2
from scipy.ndimage.filters import gaussian_filter1d
import waveform_metrics as wf_metrics



params = {
    'isi_threshold': 0.0015,
    'min_isi': 0.000166,
    'num_channels_to_compare': 7,
    'max_spikes_for_unit': 500,
    'max_spikes_for_nn': 10000,
    'n_neighbors': 4,
    'n_silhouette': 10000,
    'drift_metrics_interval_s': 51,
    'drift_metrics_min_spikes_per_interval': 10,
    'include_pc_metrics': False,
    'sample_rate': 20000,
    'maxNPeaks': 2,
    'maxNTroughs': 1,
    "minWvDuration": 0.1,  # in ms
    "maxWvDuration": 1.150,  # in ms
    "minSpatialDecaySlope": -0.008,
    "minSpatialDecaySlopeExp": 0.01,  # in a.u / um
    "maxSpatialDecaySlopeExp": 0.1,  # in a.u / um
    "maxWvBaselineFraction": 0.3,  # maximum absolute value in waveform baseline should not
    # exceed this fraction of the waveforms's absolute peak
    "maxScndPeakToTroughRatio_noise": 0.8, 
    "minTroughToPeak2Ratio_nonSomatic": 5,
    "minWidthFirstPeak_nonSomatic": 4,
    "minWidthMainTrough_nonSomatic": 5,
    "maxPeak1ToPeak2Ratio_nonSomatic": 3,
    "maxMainPeakToTroughRatio_nonSomatic": 0.8,
    ## Other classification parameters
    "minAmplitude": 40,  # in uV
    "maxRPVviolations": 0.1,  # max fraction of refractory period violations
    "maxPercSpikesMissing": 20,  # max percentage of missing spikes
    "minNumSpikes": 300,  # minimum number of total spikes recorded
    "maxDrift": 100,  # in um
    "minPresenceRatio": 0.7,  # minimum fraction of time chunks unit must be present for
    "minSNR": 5,  # min SNR for a good unit
}

FILE_PATH = r'D:\work\sorting data\2025-07-25_13-50-44-freq-300-5k-425'

def main():
    data = scio.loadmat(os.path.join(FILE_PATH,'cluster_info_per_channel.mat'))
    spike_times = data['spike_times']
    spike_clusters = data['spike_clusters']
    all_wfs = data['spikes']
    mean_wfs= data['mean_wf']

    metrics = pd.DataFrame()

    for ch in range(len(spike_clusters)):
        print(ch)

        
        ch_spike_clusters = spike_clusters[ch][0]
        ch_spike_times = spike_times[ch][0].transpose()
        ch_wfs = all_wfs[ch][0]
        ch_mean_wfs = mean_wfs[ch][0]
        duration = np.max(ch_spike_times) # 用最后一次spike时间近似
        total_units = len(np.unique(ch_spike_clusters))

        print("Calculating isi violations")
        isi_viol = calculate_isi_violations(spike_times=ch_spike_times,
                                            spike_clusters=ch_spike_clusters,
                                            total_units=total_units,
                                            isi_threshold=params['isi_threshold'],
                                            min_isi=params['min_isi'],
                                            duration=duration,
                                            verbose=True)

        print("Calculating presence ratio")
        presence_ratio = calculate_presence_ratio(spike_times=ch_spike_times,
                                                  spike_clusters=ch_spike_clusters,
                                                  total_units=total_units,
                                                  duration=duration, verbose=True)

        print("Calculating firing rate")
        firing_rate = calculate_firing_rates(ch_spike_times, ch_spike_clusters, total_units)

        print("Calculating SNR")
        snr = calculate_snrs(ch_wfs, ch_spike_clusters, total_units)

        print("Calculating amplitude cutoff")
        amplitudes = np.max(np.abs(ch_wfs),axis=0)
        amplitude_cutoff = calculate_amplitude_cutoff(spike_clusters=ch_spike_clusters,
                                                      amplitudes=amplitudes,
                                                      total_units=total_units,
                                                      verbose=True)
        
        # Waveform based metrics
        timestamps = np.arange(42) / params['sample_rate']
        wf_duration, PTratio = calculate_waveform_features(ch_mean_wfs, ch_spike_clusters, total_units, timestamps)

        metrics = pd.concat((metrics, pd.DataFrame(data=OrderedDict((('cluster_id', np.arange(total_units)),
                                                                     ('channel',ch),
                                                                     ('firing_rate', firing_rate),
                                                                     ('presence_ratio', presence_ratio),
                                                                     ('isi_violation', isi_viol),
                                                                     ('amplitude_cutoff', amplitude_cutoff),
                                                                     ('snr', snr),
                                                                     ('wf_duration', wf_duration),
                                                                     ('Peak trough ratio', PTratio),
                                                                     )))))
    return metrics
    




# ===============================================================

# HELPER FUNCTIONS TO LOOP THROUGH CLUSTERS:

# ===============================================================

def calculate_waveform_features(waveforms, spike_clusters, total_units, timestamps):
    
    durations = np.zeros((total_units,))
    PTratios = np.zeros((total_units,))

    cluster_ids = np.unique(spike_clusters)
    for idx, _ in enumerate(cluster_ids):
        durations[idx] = wf_metrics.calculate_waveform_duration(waveforms[:,idx],timestamps)
        PTratios[idx] = wf_metrics.calculate_waveform_PT_ratio(waveforms[:,idx])

    return durations,PTratios









def calculate_snrs(waveforms, spike_clusters, total_units):
    snrs = np.zeros((total_units,))
    cluster_ids = np.unique(spike_clusters)
    for idx, cluster_id in enumerate(cluster_ids):
        for_this_cluster = (spike_clusters == cluster_id)
        snrs[idx] = wf_metrics.calculate_snr(waveforms[:,for_this_cluster.squeeze()])

    return snrs




def calculate_isi_violations(spike_times, spike_clusters, total_units, isi_threshold, min_isi, duration,
                             spike_cluster_subset=None, verbose=True):
    if spike_cluster_subset is not None:
        cluster_ids = spike_cluster_subset
    else:
        cluster_ids = np.unique(spike_clusters)

    viol_rates = np.zeros((total_units,))

    for idx, cluster_id in enumerate(cluster_ids):


        for_this_cluster = (spike_clusters == cluster_id)
        viol_rates[cluster_id-1], num_violations = isi_violations(spike_times[for_this_cluster],
                                                                duration=duration,
                                                                isi_threshold=isi_threshold,
                                                                min_isi=min_isi)

    return viol_rates


def calculate_presence_ratio(spike_times, spike_clusters, total_units, duration, spike_cluster_subset=None,
                             verbose=True):
    if spike_cluster_subset is not None:
        cluster_ids = spike_cluster_subset
    else:
        cluster_ids = np.unique(spike_clusters)

    ratios = np.zeros((total_units,))

    for idx, cluster_id in enumerate(cluster_ids):


        for_this_cluster = (spike_clusters == cluster_id)
        ratios[cluster_id-1] = presence_ratio(spike_times[for_this_cluster],
                                            duration=duration)

    return ratios


def calculate_num_spikes(spike_times, spike_clusters, total_units, spike_cluster_subset=None, verbose=True):
    num_spikes = np.zeros((total_units,))
    if spike_cluster_subset is not None:
        cluster_ids = spike_cluster_subset
    else:
        cluster_ids = np.unique(spike_clusters)

    for idx, cluster_id in enumerate(cluster_ids):


        for_this_cluster = (spike_clusters == cluster_id)
        num_spikes[cluster_id-1] = len(spike_times[for_this_cluster])

    return num_spikes


def calculate_firing_rates(spike_times, spike_clusters, total_units):

    cluster_ids = np.unique(spike_clusters)

    firing_rates = np.zeros((total_units,))

    min_time = np.min(spike_times)
    max_time = np.max(spike_times)

    for idx, cluster_id in enumerate(cluster_ids):


        for_this_cluster = (spike_clusters == cluster_id)
        firing_rates[idx] = firing_rate(spike_times[for_this_cluster],
                                        min_time = np.min(spike_times),
                                        max_time = np.max(spike_times))

    return firing_rates



def calculate_amplitude_cutoff(spike_clusters, amplitudes, total_units, spike_cluster_subset=None, verbose=True):
    if spike_cluster_subset is not None:
        cluster_ids = spike_cluster_subset
    else:
        cluster_ids = np.unique(spike_clusters)

    amplitude_cutoffs = np.zeros((total_units,))

    for idx, cluster_id in enumerate(cluster_ids):


        for_this_cluster = (spike_clusters == cluster_id)
        amplitude_cutoffs[cluster_id-1] = amplitude_cutoff(amplitudes[for_this_cluster.squeeze()])

    return amplitude_cutoffs


# ==========================================================

# IMPLEMENTATION OF ACTUAL METRICS:

# ==========================================================

def firing_rate(spike_train, min_time = None, max_time = None):
    """Calculate firing rate for a spike train.

    If no temporal bounds are specified, the first and last spike time are used.

    Inputs:
    -------
    spike_train : numpy.ndarray
        Array of spike times in seconds
    min_time : float
        Time of first possible spike (optional)
    max_time : float
        Time of last possible spike (optional)

    Outputs:
    --------
    fr : float
        Firing rate in Hz

    """

    if min_time is not None and max_time is not None:
        duration = max_time - min_time
    else:
        duration = np.max(spike_train) - np.min(spike_train)

    fr = spike_train.size / duration

    return fr


def isi_violations(spike_train, duration, isi_threshold, min_isi=0):
    """Calculate Inter-Spike Interval (ISI) violations for a spike train.

    Based on metric described in Hill et al. (2011) J Neurosci 31: 8699-8705

    Originally written in Matlab by Nick Steinmetz (https://github.com/cortex-lab/sortingQuality)
    Converted to Python by Daniel Denman

    Inputs:
    -------
    spike_train : array of monotonically increasing spike times (in seconds) [t1, t2, t3, ...]
    duration : length of recording (seconds)
    isi_threshold : threshold for classifying adjacent spikes as an ISI violation
      - this is the biophysical refractory period
    min_isi : minimum possible inter-spike interval (default = 0)
      - this is the artificial refractory period enforced by the data acquisition system
        or post-processing algorithms

    Outputs:
    --------
    fpRate : rate of contaminating spikes as a fraction of overall rate
      - higher values indicate more contamination
    num_violations : total number of violations detected

    """
    isis_initial = np.diff(spike_train)

    if min_isi > 0:
        duplicate_spikes = np.where(isis_initial <= min_isi)[0]
        spike_train = np.delete(spike_train, duplicate_spikes + 1)

    isis = np.diff(spike_train)
    num_spikes = len(spike_train)
    num_violations = sum(isis < isi_threshold)
    violation_time = 2 * num_spikes * (isi_threshold - min_isi)
    total_rate = firing_rate(spike_train, duration)
    violation_rate = num_violations / violation_time
    fpRate = violation_rate / total_rate

    return fpRate, num_violations


def presence_ratio(spike_train, duration, num_bin_edges=101):
    """Calculate fraction of time the unit is present within an epoch.

    Inputs:
    -------
    spike_train : array of spike times
    duration : length of recording (seconds)
    num_bin_edges : number of bin edges for histogram
      - total bins = num_bin_edges - 1

    Outputs:
    --------
    presence_ratio : fraction of time bins in which this unit is spiking

    """

    h, b = np.histogram(spike_train, np.linspace(0, duration, num_bin_edges))

    return np.sum(h > 0) / (num_bin_edges - 1)



def amplitude_cutoff(amplitudes, num_histogram_bins=500, histogram_smoothing_value=3):
    """ Calculate approximate fraction of spikes missing from a distribution of amplitudes

    Assumes the amplitude histogram is symmetric (not valid in the presence of drift)

    Inspired by metric described in Hill et al. (2011) J Neurosci 31: 8699-8705

    Input:
    ------
    amplitudes : numpy.ndarray
        Array of amplitudes (don't need to be in physical units)
    num_histogram_bins : int
        Number of bins for calculating amplitude histogram
    histogram_smoothing_value : float
        Gaussian filter window for smoothing amplitude histogram

    Output:
    -------
    fraction_missing : float
        Fraction of missing spikes (ranges between 0 and 0.5)
        If more than 50% of spikes are missing, an accurate estimate isn't possible

    """

    h, b = np.histogram(amplitudes, num_histogram_bins, density=True)

    pdf = gaussian_filter1d(h, histogram_smoothing_value)
    support = b[:-1]

    peak_index = np.argmax(pdf)
    G = np.argmin(np.abs(pdf[peak_index:] - pdf[0])) + peak_index

    bin_size = np.mean(np.diff(support))
    fraction_missing = np.sum(pdf[G:]) * bin_size

    fraction_missing = np.min([fraction_missing, 0.5])

    return fraction_missing


def make_index_mask(spike_clusters, unit_id, min_num, max_num, seed=None):
    """ Create a mask for the spike index dimensions of the pc_features array

    Inputs:
    -------
    spike_clusters : numpy.ndarray (num_spikes x 0)
        Contains cluster IDs for all spikes in pc_features array
    unit_id : Int
        ID for this unit
    min_num : Int
        Minimum number of spikes to return; if there are not enough spikes for this unit, return all False
    max_num : Int
        Maximum number of spikes to return; if too many spikes for this unit, return a random subsample
    seed: int
        Random seed for reproducibility

    Output:
    -------
    index_mask : numpy.ndarray (boolean)
        Mask of spike indices for pc_features array

    """

    index_mask = spike_clusters == unit_id

    inds = np.where(index_mask)[0]

    if len(inds) < min_num:
        index_mask = np.zeros((spike_clusters.size,), dtype='bool')
    else:
        index_mask = np.zeros((spike_clusters.size,), dtype='bool')
        order = np.random.RandomState(seed=seed).permutation(inds.size)
        index_mask[inds[order[:max_num]]] = True

    return index_mask

if __name__ == '__main__':
    metrics = main()
    # save to csv
    metrics.to_csv('metrics.csv', index=False)


