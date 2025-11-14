#!/usr/bin/env python
# coding: utf-8

# In[16]:


import sys
import os

# Añade el directorio padre al sys.path
sys.path.append(os.path.abspath("/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/WCTE_BRB_Data_Analysis"))
sys.path.append(os.path.abspath("/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/nHits_trigger"))
sys.path.append('/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/hipy')
sys.path.append("/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/WCTE_event_display")

import hipy.pltext       as pltext
import hipy.utils        as ut
import matplotlib.pyplot as plt
import awkward           as ak
import numpy             as np
import matplotlib.colors as colors

from src.read_data import load_concatenated, read_parquet, nHits
from tqdm          import tqdm
from EventDisplay  import EventDisplay

pltext.style()


# ## Read Data

# In[17]:


# ============= LOADS THE DATA ===================
bkg_run  = 1766
bkg_data = load_concatenated(outdir=f"/mnt/lustre/scratch/nlsas/home/usc/ie/dcr/hk/processed_wcte_data/{bkg_run}_1_calibrated_tof")
print("Data Loaded")

# =============== MASK OUT CARDS 130-132 ===================
bkg_run_cards, bkg_run_channels, bkg_run_slots, bkg_run_positions, bkg_run_times, bkg_run_events, bkg_run_charges, bkg_run_window_times = read_parquet(bkg_data, mask=True)
print("Data Masked")

bkg_run_times_sorted = ak.sort(bkg_run_times) - bkg_run_window_times
bkg_run_window_times_sorted = ak.sort(bkg_run_window_times, axis=0)


# In[18]:


# ============= LOADS THE DATA ===================
nicf_run  = 1767
nicf_data = load_concatenated(outdir=f"/mnt/lustre/scratch/nlsas/home/usc/ie/dcr/hk/processed_wcte_data/{nicf_run}_1_calibrated_tof")
print("Data Loaded")

# =============== MASK OUT CARDS 130-132 ===================
nicf_run_cards, nicf_run_channels, nicf_run_slots, nicf_run_positions, nicf_run_times, nicf_run_events, nicf_run_charges, nicf_run_window_times = read_parquet(nicf_data, mask=True)
print("Data Masked")

nicf_run_times_sorted = ak.sort(nicf_run_times) - nicf_run_window_times
nicf_run_window_times_sorted = ak.sort(nicf_run_window_times, axis=0)


# # nHits To Find Spills and AfterPulsing

# In[19]:


bkg__triggered_spill_hits_index, bkg__triggered_spill_hit_times = nHits(mode="multiple_events", hit_times=bkg_run_times_sorted, w=5000, thresh_min=300, thresh_max=10000, pre_window=0, post_window=4000, jump=9000)
nicf_triggered_spill_hits_index, nicf_triggered_spill_hit_times = nHits(mode="multiple_events", hit_times=nicf_run_times_sorted, w=5000, thresh_min=300, thresh_max=10000, pre_window=0, post_window=4000, jump=9000)


# ### Remove Spills from data

# In[20]:


def remove_spill(trigger_indices, run_times):
    corrected_run_times = []

    for i in tqdm(range(len(run_times))):
        if i in trigger_indices:
            #Event with spill → remove selected indices
            data = run_times[i]
            trigger_indices_event = np.concatenate(trigger_indices[i])
            all_indices           = np.arange(len(data))
            valid_indices         = np.setdiff1d(all_indices, trigger_indices_event)
            corrected_run_times.append(data[valid_indices])
        else:
            # Event without spill → leave as it is
            corrected_run_times.append(run_times[i])

    return ak.Array(corrected_run_times)


# In[21]:


bkg__corrected_run_times_sorted = remove_spill(bkg__triggered_spill_hits_index, bkg_run_times_sorted)
nicf_corrected_run_times_sorted = remove_spill(nicf_triggered_spill_hits_index, nicf_run_times_sorted)


# In[ ]:


# bkg__corrected_no_hits_per_bin_50ns = np.concatenate([np.histogram(bkg__corrected_run_times_sorted[ev], np.arange(0, 500e3, 50))[0] for ev in tqdm(range(len(bkg__corrected_run_times_sorted)), total=len(bkg__corrected_run_times_sorted))])
# nicf_corrected_no_hits_per_bin_50ns = np.concatenate([np.histogram(nicf_corrected_run_times_sorted[ev], np.arange(0, 500e3, 50))[0] for ev in tqdm(range(len(nicf_corrected_run_times_sorted)), total=len(nicf_corrected_run_times_sorted))])

bkg__corrected_no_hits_per_bin_50ns = np.concatenate([np.histogram(bkg__corrected_run_times_sorted[ev], np.arange(0, 500e3, 50))[0] for ev in tqdm(range(1000), total=1000)])
nicf_corrected_no_hits_per_bin_50ns = np.concatenate([np.histogram(nicf_corrected_run_times_sorted[ev], np.arange(0, 500e3, 50))[0] for ev in tqdm(range(1000), total=1000)])

bkg__corrected_no_hits_per_bin_50ns = [i for i in bkg__corrected_no_hits_per_bin_50ns if i != 0]
nicf_corrected_no_hits_per_bin_50ns = [i for i in nicf_corrected_no_hits_per_bin_50ns if i != 0]

# In[ ]:


subplot = pltext.canvas(2)
bin_width = 50

subplot(1)
c501, b501, _501 = pltext.hist(bkg__corrected_no_hits_per_bin_50ns, 100, ylog=True, range=(0, 300), density=True, xylabels=f"No. Hits Per {bin_width} ns bin", label="bkg");
c502, b502, _502 = pltext.hist(nicf_corrected_no_hits_per_bin_50ns, 100, ylog=True, range=(0, 300), density=True, xylabels=f"No. Hits Per {bin_width} ns bin", label="nicf");
plt.savefig("1.pdf")

subplot(2)
pltext.hist(bkg__corrected_no_hits_per_bin_50ns, 19, ylog=True, range=(2, 20), xylabels=f"No. Hits Per {bin_width} ns bin", label="bkg");
pltext.hist(nicf_corrected_no_hits_per_bin_50ns, 19, ylog=True, range=(2, 20), xylabels=f"No. Hits Per {bin_width} ns bin", label="nicf");
plt.savefig("2.pdf")


# In[ ]:




