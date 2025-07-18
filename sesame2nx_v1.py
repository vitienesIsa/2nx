#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 12:37:24 2025

@author: adminisa
"""

# VERSION 1

# Developed in 06.2025 by Isabela Vitienes, in the Zaslansky lab 
# of the Department of Operative, Preventive and Pediatric Dentistry 
# at Charité – Universitätsmedizin Berlin, with funding from the 
# DFG (FOR5657).


##################################################################
# This code converts tomography datasets collected at Sesame     #
# into NeXus (.nx) format which is compatible with Nabu tools    # 
# for image reconstruction.                                      #
##################################################################

#########################################################
#                                                       #  
# INPUT: (1) User-selected .h5 file from Sesame BEATS   #
#        (2) User-selected output directory             #
#        (3) Pixel size, in um                          #
#        (4) Beam energy, in keV                        #
#        (5) Binning factor                             #
# OUTPUT: Version of input data in NeXus (.nx) format   #
#                                                       #
#########################################################



import numpy as np
import h5py
from nxtomo import NXtomo
from nxtomo.nxobject.nxdetector import ImageKey
import os
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from datetime import datetime
import time
from tqdm import tqdm
import sys
import threading
import itertools




###############################################################################
# functions ###################################################################
###############################################################################

def spinner_task(stop_event):
    spinner = itertools.cycle(['|', '/', '-', '\\'])
    while not stop_event.is_set():
        sys.stdout.write(next(spinner))
        sys.stdout.flush()
        time.sleep(0.1)
        sys.stdout.write('\b')

def upstream_path (path):
    for i in range(len(path)-2, -1, -1):
        if path[i] == '/':
            k = i
            break
    return path[:k+1], path[k+1:]



###############################################################################
# user inputs #################################################################
###############################################################################

# Input file path
root = tk.Tk()
root.withdraw()  # Hide the root window
input_path = filedialog.askopenfilename(title="Select an .h5 file", filetypes=[("HDF5 files", "*.h5")])
filename = upstream_path(input_path)[1][:-3]
print('You are now converting the file ' + input_path)
fbam = h5py.File(input_path, 'r') # Read as h5 object
z, x, y = fbam['exchange/data'].shape
print('Your dataset has ' + str(z-1) + ' projections.')
print('Each projection is ' + str(x) + ' x ' + str(y) + ' pixels.')

scan_type = 'full'

# Prompt user to select output directory
root = tk.Tk()
root.withdraw()  # Hide the root window
op = filedialog.askdirectory(title="Select output location in /data")
print('The output .nx file will be stored here: ' + op)


# Prompt user to input pixel size
# Seems to be some error in the sesame raw data files, where pixel size always says 6.5
root = tk.Tk()
root.withdraw()
pixel_size = simpledialog.askfloat(title="Pixel size", prompt="Specify pixel size in units of um")
print('Pixel size: ' + str(pixel_size) + ' um')

# Prompt user to input beam energy
root = tk.Tk()
root.withdraw()
energy = simpledialog.askinteger(title="Energy", prompt="Specify beam energy in units of keV")
energy = float(energy)
print('Beam energy: ' + str(energy) + ' keV')


# Prompt user to input binning factor, i.e. number of projections to average. 
# Ideally will be a factor of the total number of projections. Otherwise, the last averaged set will be smaller than the rest. 
root = tk.Tk()
root.withdraw()
proj_bin_size = simpledialog.askinteger(title="Projection binning", prompt="How many projections do you want to average?")
print('Every ' + str(proj_bin_size) + ' projections will be averaged.')


###############################################################################
# computation #################################################################
###############################################################################

start_time = datetime.now()
print('FYI: Conversion may take a while. How long depends on the size of your data and available resources. ')
print('Start = ' + start_time.strftime('%Hh:%Mm:%Ss, %m.%d.%Y'))



# Extract dataset from path
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Loading input file... ", end="", flush=True)

stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Input file loaded')


# Initialize nx output
new_nx = NXtomo() 


# Extract all projections
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Extracting projections... ", end="", flush=True)
data = np.array(fbam['exchange/data'])[1:]
n = np.shape(data)[0]
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Projections extracted')
print('Number of projections: ' + str(len(data)))


# Extract and reformat flats
# Note: .nx has has placeholding for two sets of flats. If only one taken, a copy of this one used in place of second set.
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Extracting and reformatting flats and darks... ", end="", flush=True)
# =====================================
flats1_stack = np.array(fbam['exchange/data_white'])
flats2_stack = np.array(fbam['exchange/data_white'])
darks_stack = np.array(fbam['exchange/data_dark'])
# =====================================
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Flats and darks extracted and reformatted')


# Extract 'true' projections and simultaneously bin
n_truradio = len(data)
if proj_bin_size == 1:
    no_bins = n_truradio
else:
    no_bins = int(n_truradio/proj_bin_size) + 1 


if proj_bin_size == 1:
    stop_spinner = threading.Event()
    spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
    spinner_thread.start()
    print("Extracting projections. May take a while, dont despair.", end="", flush=True)
    binned_projs = data
    stop_spinner.set()
    spinner_thread.join()
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Projections extracted.')
else:
    for i in tqdm(range(no_bins), desc='Extracting and binning projections. May take some time, dont despair.'):
        a = i * proj_bin_size
        b = (i+1)*proj_bin_size - 1
        if i == 0: # first bin, initialize binned data variable      
            binned_projs = np.mean(data[a : b], axis=0)
            binned_projs = binned_projs.reshape(1, binned_projs.shape[0], binned_projs.shape[1])
        elif i != no_bins-1: # not last bin
            avg = np.mean(data[a : b], axis=0)
            binned_projs = np.vstack([binned_projs, [avg]])
        elif i == no_bins-1: # last bin
            avg = np.mean(data[a :], axis=0)
            binned_projs =  np.vstack([binned_projs, [avg]])
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Projections extracted and binned') 
print('Number of binned projections: ' + str(len(binned_projs)))



# Reformatting projections
radiosA = binned_projs[:int(no_bins/2)]
radiosB = binned_projs[int(no_bins/2):]
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Binned projections reformatted')
print('first half of radios: ' + str(len(radiosA)))
print('second half of radios: ' + str(len(radiosB)))


# BAM doesnt take alignment images, so we just take the halfway projection, i.e. first radiograph in radiosB
alignment_stack = radiosB[0,:,:]
alignment_stack = np.reshape(alignment_stack, (1, alignment_stack.shape[0], alignment_stack.shape[1]))


# Ensure all stacks are of right dimension
assert darks_stack.ndim == 3
assert flats1_stack.ndim == 3
assert flats2_stack.ndim == 3
assert radiosA.ndim == 3
assert radiosB.ndim == 3
assert alignment_stack.ndim == 3



# Compile all data, flats darks radios alignment, into a dataset, in the order specified bx NeXus convention
new_data = np.concatenate([darks_stack, flats1_stack, radiosA, flats1_stack, radiosB, flats2_stack, alignment_stack])
print(len(new_data))
assert new_data.ndim == 3


# Add this new (reorganized) dataset to the new nx object
new_nx.instrument.detector.data = new_data
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Reformatted data compiled into nx object')


# Image key control
# Metadata to give identity to each slice in the data
image_key_control = np.concatenate([
    [ImageKey.DARK_FIELD] * int(len(darks_stack)),
    [ImageKey.FLAT_FIELD] * int(len(flats1_stack)),
    [ImageKey.PROJECTION] * int(len(radiosA)),
    [ImageKey.FLAT_FIELD] * int(len(flats1_stack)),
    [ImageKey.PROJECTION] * int(len(radiosB)),
    [ImageKey.FLAT_FIELD] * int(len(flats2_stack)),
    [ImageKey.ALIGNMENT] * int(len(alignment_stack)),
])
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Image key control created')
print('len image key control: ' + str(len(image_key_control)))
assert len(image_key_control) == len(new_data)
new_nx.instrument.detector.image_key_control = image_key_control



# Rotation angles
# Metadata specifying rotation angle for each slice
angles_raw = np.array(fbam['exchange/theta'])[1:]
scanspan = int(np.max(angles_raw))
angles = np.linspace(0,scanspan,no_bins)
rotation_angle = np.concatenate([
    [0.0, ] * int(len(darks_stack)),
    [0.0, ] * int(len(flats1_stack)),
    angles[:int(len(angles)/2)],
    [scanspan/2, ] * int(len(flats1_stack)),
    angles[int(len(angles)/2):],
    [scanspan, ] * int(len(flats2_stack)),
    [scanspan/2, ] * int(len(alignment_stack)),
])
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Rotation angles specified')
print('len of rotation angles: ' + str(len(rotation_angle)))

assert len(rotation_angle) == len(new_data)
new_nx.sample.rotation_angle = rotation_angle



# Field of view
new_nx.instrument.detector.field_of_view = str(scan_type) 

# Pixel size
new_nx.instrument.detector.x_pixel_size = new_nx.instrument.detector.y_pixel_size = pixel_size*1e-6  # pixel size in SI: meter

# Energy and detector distance, needed for phase retrieval
new_nx.energy = energy

detector_distance = fbam['measurement/instrument/detector_motor_stack/detector_z'][0] # in units of mm
new_nx.instrument.detector.distance = detector_distance*1e-4 # in units of m


ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Additional scan metadata specified (pixel size, detector distance, etc.)')


root_out_path = op + '/' + filename + '/'
if not os.path.exists(root_out_path):
    os.mkdir(root_out_path)
out_path = root_out_path + filename[:-3] + '_bin' + str(proj_bin_size) + '_v2_nxtomo.nx'
if os.path.exists(out_path):
    out_path = out_path[:-3] + '_1.nx'
nx_tomo_file_path = os.path.join (out_path)
new_nx.save(file_path=nx_tomo_file_path, data_path='entry', overwrite=True)
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Out nx file saved')

end_time = datetime.now()
print('End = ' + end_time.strftime('%Hh:%Mm:%Ss, %m.%d.%Y'))
duration = start_time - end_time
total_seconds = int(duration.total_seconds())
hours, remainder = divmod(total_seconds, 3600)
minutes, seconds = divmod(remainder, 60)
print(f'Nx conversion COMPLETE. Conversion duration = {hours}h:{minutes}m:{seconds}s')
print('Find your converted file in ' + out_path)

