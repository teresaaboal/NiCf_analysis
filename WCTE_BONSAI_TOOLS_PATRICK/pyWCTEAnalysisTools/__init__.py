import ROOT 
import glob 
import numpy as np 
import pandas as pd 
import array
import pandas as pd
import os
import cppyy

def get_offline_run_base(run, base="/eos/experiment/wcte/data/2025_commissioning/"):
    return f"{base}/processed_offline_data/production_v0/{run}/WCTE_offline_R*.root"

def get_offline_run_files(run, base="/eos/experiment/wcte/data/2025_commissioning/"):
    files = glob.glob(f"{base}/processed_offline_data/production_v0/{run}/WCTE_offline_R*.root")
    return files

def get_offline_run_tchain(run, limit=0):
    files = get_offline_run_files(run)
    tt = ROOT.TChain("WCTEReadoutWindows")
    count = 0
    for ff in files:
        tt.Add(ff)
        count += 1
        if limit > 0 and count > limit: break
    
    return tt

# Function to load geometry mapping from a text file
def get_geo_mapping():
    geo = pd.read_csv("geofile_NuPRISMBeamTest_16cShort_mPMT.txt", 
                      index_col=False,      # Do not use any column as the index
                      sep='\s+',            # Use whitespace as separator
                      skiprows=5,           # Skip the first 5 header lines
                      names=["id","mpmtid","spmtid",
                             "x","y","z","dx","dy","dz", "cyloc"])  # Explicit column names

    return geo


# Global Geo definition : Change to MC
geo = get_geo_mapping()


def getxyz(geo, mpmtids, posids):
    # Build a single lookup dictionary: {(mpmtid, spmtid): (x, y, z, id)}
    lookup = {
        (row.mpmtid, row.spmtid): (row.x, row.y, row.z, row.id)
        for row in geo.itertuples(index=False)
    }

    # Adjust input IDs to match geometry convention
    keys = zip((mid + 1 for mid in mpmtids), (sid + 1 for sid in posids))

    # Use the lookup dictionary to retrieve values efficiently
    results = [lookup.get(k, (-999.9, -999.9, -999.9, -999)) for k in keys]

    # Unpack results into separate arrays
    if len(results) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    
    x, y, z, c = map(np.array, zip(*results))

    return x, y, z, c



# Simple class to store a filtered time slice of hit data
class hit_slice:
    def __init__(self, hits):
        self.t = hits.t  # Hit times
        self.q = hits.q  # Hit charges
        self.x = hits.x  # x-positions
        self.y = hits.y  # y-positions
        self.z = hits.z  # z-positions
        self.cable = hits.cable 

# Class to process and store a collection of hits from an event
class hit_collection:
    def __init__(self, e, g=geo):
        
        # Apply time window cut for hit calibration times
        above_10us = True #np.array(e.hit_pmt_calibrated_times) >   10000
        below_490us = True # np.array(e.hit_pmt_calibrated_times) < 490000
        has_calib = np.array(e.hit_pmt_has_time_constant)     # Check for valid calibration
        
        has_valid_card = ((np.array(e.hit_mpmt_slot_ids) <= 106) & 
                          (np.array(e.hit_pmt_position_ids) <= 19) & 
                          (np.array(e.hit_mpmt_card_ids) < 120))
        
        has_valid_charge = ((np.array(e.hit_pmt_charges) < 1e4))
                            
        hit_selection = has_valid_charge & above_10us & below_490us & has_calib & has_valid_card # Combined mask
        
        # Apply selection mask to extract relevant hit information
        self.t = np.array(e.hit_pmt_calibrated_times)[hit_selection]
        self.q = np.array(e.hit_pmt_charges)[hit_selection]
        
        self.card = np.array(e.hit_mpmt_card_ids)[hit_selection]
        self.slot = np.array(e.hit_mpmt_slot_ids)[hit_selection]
        self.channel = np.array(e.hit_pmt_channel_ids)[hit_selection]
        self.posid = np.array(e.hit_pmt_position_ids)[hit_selection]
        
        # Get spatial positions and cable IDs using geometry mapping
        self.x, self.y, self.z, self.cable = getxyz(g, self.slot, self.posid)
        
        # Remove hits with invalid geometry mapping (sentinel x = -999.9)
        cut=self.x > -999
        self.t = self.t[cut]
        self.q = self.q[cut]
        self.card = self.card[cut]
        self.slot = self.slot[cut]
        self.channel = self.channel[cut]
        self.posid = self.posid[cut]
        self.z = self.z[cut]
        self.y = self.y[cut]
        self.cable = self.cable[cut]
        self.x = self.x[cut]
        
    # Return a new hit_slice object within a time window
    def time_slice(self, tstart, tend):
        h = hit_slice(self)  # Create a new slice from this collection
        
        tsel = (h.cable > -1) & (h.t > tstart) & (h.t < tend)  # Select hits within the time range
        
        # Apply selection to all relevant hit attributes
        h.q = h.q[ tsel ]
        h.x = h.x[ tsel ]
        h.y = h.y[ tsel ]
        h.z = h.z[ tsel ]
        h.cable = h.cable[ tsel ]
        h.t = h.t[ tsel ]
                
        return h  # Return the filtered slice


    
    
    
    
    
    