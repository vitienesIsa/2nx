#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 12:37:24 2025

@author: adminisa
"""

# VERSION 1

# Developed in the Zaslansky lab of the Department of Operative, 
# Preventive and Pediatric Dentistry at Charité – Universitätsmedizin 
# Berlin, with funding from the DFG (FOR5657).

##################################################################
# This code converts tomography datasets collected at BAMline    #
# into NeXus (.nx) format which is compatible with Nabu tools    # 
# for image reconstruction.                                      #
##################################################################

#########################################################
#                                                       #  
# INPUT: (1) User-selected .h5 file from BESSY BAMline  #
#        (2) User-selected output directory             #
#        (3) Binning factor                             #
# OUTPUT: Version of input data in NeXus (.nx) format   #
#                                                       #
#########################################################



import numpy as np
import h5py
from nxtomo import NXtomo
from nxtomo.nxobject.nxdetector import ImageKey
import statistics as stat
import os
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import ttk
from datetime import datetime
from tqdm import trange
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
z, x, y = fbam['entry/data/data'].shape
print('Your dataset has ' + str(z-1) + ' projections.')
print('Each projection is ' + str(x) + ' x ' + str(y) + ' pixels.')

scan_type = 'full'

# Prompt user to select output directory
root = tk.Tk()
root.withdraw()  # Hide the root window
op = filedialog.askdirectory(title="Select output location in /data")
print('The output .nx file will be stored here: ' + op)

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


# Initialize nx output
new_nx = NXtomo() 


# Extract all projections
data = fbam['entry/data/data'] 
n = np.shape(data)[0]
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Input data loaded')
print('Length of full data: ' + str(len(data)))


# timestamp metadata
ndarraytimestamp = np.array(fbam['entry/instrument/NDAttributes/NDArrayTimeStamp']) # take advantage of fact that there is a larger increase in timestep at end/beginning of flat acquisition
dd_ndarraytimestep = np.array(np.diff(ndarraytimestamp)) # timestamp increments 
mode = stat.mode(dd_ndarraytimestep) # most common timestamp increment, assume to be increment once a stable velocity is achieved


# Find locations (in terms of indices of data variable) of flats
k = []
for i in range(len(dd_ndarraytimestep)):
    if dd_ndarraytimestep[i] > mode + 3:
        k = k + [i]
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Flats located')
print('Index locations of flats: ' + str(k))


# Extract and reformat flats
# Note: .nx has has placeholding for two sets of flats. If only one taken, a copy of this one used in place of second set.
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Extracting and reformatting flats... ", end="", flush=True)
# =====================================
if len(k) == 1: # only one set of flats
    no_flats = [[k[0] + 1],[]]
    if k[0] < 100: # flats taken before radios (flat1_stack). assume there will never be more than 100 flats, or less than 200 radios
        f1_end_ind = k[0] + 1
        flats1_stack = np.array(data[:f1_end_ind])
        flats2_stack = flats1_stack # since no flats taken after radios
        k = k[0] + [0] # adjust k to reflect that there are zero flats taken after radios
    elif k[0] > 100: # flats taken after radios (flat2_stats)
        f2_start_ind = k[0] + 1
        flats2_stack = np.array(data[f2_start_ind:])
        flats1_stack = flats2_stack # since no flats taken before radios
        k = [0] + k[0] # adjust k to reflect that there are zero flats taken before radios
if len(k) == 2: # two sets of flats, one before and one after radios
    no_flats = [[k[0] + 1], [n - k[1] - 1]]
    f1_end_ind = k[0] + 1
    f2_start_ind = k[1] + 1
    flats1_stack = np.array(data[:f1_end_ind])
    flats2_stack = np.array(data[f2_start_ind:])
# =====================================
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()

print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Flats extracted and reformatted')


# BAM doesnt produce darks so we just use 'zero' arrays
darks_stack = np.array(np.zeros_like(flats1_stack))+100 
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Darks created')


# Identify and remove accelerating radios
# In bam, motor initially accelerates until reaching target speed, and radios taken during this acceleration must be removed 
# Look to 'SAMPLE_W' data, which gives the rotation angle. 
samplew = np.array(fbam['/entry/instrument/NDAttributes/SAMPLE_W']) # rotation angles, all
samplew_radio = samplew[len(flats1_stack):-len(flats2_stack)]
scanspan = int(10*(np.floor(np.max(samplew_radio)/10))) # whether 180 or 360 scan
accspan = abs((np.max(samplew) - np.min(samplew))-scanspan) # max scan angle of acceleration tomos
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Acceleration radiographs identified')



n_allradio = n - no_flats[0][0] - no_flats[1][0] # Number of projections (number of total projections minus flats, since bam doesnt take darks)
n_acc = n_allradio%100 
n_truradio = n_allradio - n_acc # number of 'true', non-accelerating projections. Should be a 'nice' number i.e. multiple of 100.
print('Number of all projections (all radiographs excluding flats and darks): ' + str(n_allradio))
print('Number of final projections (after removing accelleration projections): ' + str(n_truradio))
shift2tru = no_flats[0][0]+n_acc

# Extract 'true' projections and simultaneously bin
if n_truradio%proj_bin_size == 0:
    no_bins = n_truradio
else:
    no_bins = int(n_truradio/proj_bin_size) + 1 

if proj_bin_size == 1:
    stop_spinner = threading.Event()
    spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
    spinner_thread.start()
    print("Extracting projections. May take a while...dont despair.", end="", flush=True)
    binned_projs = data[shift2tru:shift2tru+n_truradio]
    stop_spinner.set()
    spinner_thread.join()
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Projections extracted')
else:
    for i in tqdm(range(no_bins), desc='Extracting and binning projections. May take a while...dont despair.'):
        a = i * proj_bin_size
        b = (i+1)*proj_bin_size - 1
        if i == 0: # first bin, initialize binned data variable      
            binned_projs = np.mean(data[shift2tru + a : shift2tru + b], axis=0)
            binned_projs = binned_projs.reshape(1, binned_projs.shape[0], binned_projs.shape[1])
        elif i != no_bins-1: # not last bin
            avg = np.mean(data[shift2tru + a : shift2tru + b], axis=0)
            binned_projs = np.vstack([binned_projs, [avg]])
        elif i == no_bins-1: # last bin
            avg = np.mean(data[shift2tru + a :], axis=0)
            binned_projs =  np.vstack([binned_projs, [avg]])
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Projections extracted and binned') 
print('Number of projections: ' + str(len(binned_projs)))



# Reformatting projections
radiosA = binned_projs[:int(no_bins/2)]
radiosB = binned_projs[int(no_bins/2):]
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Projections reformatted')
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
new_nx.instrument.detector.field_of_view = str(scan_type) # Need to find out how to figure out from bam file what fov was

# Pixel size
pixel_size = fbam['entry/instrument/NDAttributes/CT_Pixelsize'][()][0] # in units of microns
new_nx.instrument.detector.x_pixel_size = new_nx.instrument.detector.y_pixel_size = pixel_size*1e-6  # pixel size in SI: meter

# Energy and detector distance, needed for phase retrieval
energy = fbam['entry/instrument/NDAttributes/DMM_Energy'][()][0]
new_nx.energy = energy

detector_distance = (fbam['entry/instrument/NDAttributes/CT-Kamera-Z'][()][0]+25) # in units of mm
new_nx.instrument.detector.distance = detector_distance*1e-4 # in units of m
#new_nx.instrument.sample_detector_distance = detector_distance*1e-4 # in units of m

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

