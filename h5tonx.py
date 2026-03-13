##################################################################
# This code converts tomography datasets collected at BAMline    #
# or Sesame BEATS into NeXus (.nx) format compatible with Nabu  #
# tools for image reconstruction.                               #
#                                                                #
# MEMORY-EFFICIENT: Data is written chunk-by-chunk so the full   #
# dataset is never held in RAM.                                  #
##################################################################

#########################################################
#                                                       #
# INPUT: (1) User-selected .h5 file                     #
#        (2) Source beamline (BAM or Sesame)            #
#        (3) User-selected output directory             #
#        (4) Binning factor                             #
#        (5) Scan type (full or half)                   #
#        (6) [Sesame only] Pixel size, energy           #
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
from tkinter import filedialog, simpledialog, messagebox
from datetime import datetime
from tqdm import tqdm
import pint
from importlib.metadata import version as pkg_version


###############################################################################
# functions ###################################################################
###############################################################################

def upstream_path(path):
    for i in range(len(path) - 2, -1, -1):
        if path[i] == '/':
            k = i
            break
    return path[:k + 1], path[k + 1:]


###############################################################################
# user inputs #################################################################
###############################################################################

# --- Input file ---
root = tk.Tk()
root.withdraw()
input_path = filedialog.askopenfilename(title="Select an .h5 file", filetypes=[("HDF5 files", "*.h5")])
filename = upstream_path(input_path)[1][:-3]
print('You are now converting the file ' + input_path)
fbam = h5py.File(input_path, 'r')

# --- Source beamline ---
root = tk.Tk()
root.withdraw()
source = simpledialog.askstring(title="Source beamline", prompt="Which beamline? Enter 'BAM' or 'Sesame'")
source = source.strip().upper()
assert source in ('BAM', 'SESAME'), "Source must be 'BAM' or 'Sesame'"
print('Source beamline: ' + source)

# --- Read dataset shape ---
if source == 'BAM':
    z, x, y = fbam['entry/data/data'].shape
else:
    z, x, y = fbam['exchange/data'].shape
print('Your dataset has ' + str(z - 1) + ' projections.')
print('Each projection is ' + str(x) + ' x ' + str(y) + ' pixels.')

# --- Output directory ---
root = tk.Tk()
root.withdraw()
op = filedialog.askdirectory(title="Select output location")
print('The output .nx file will be stored here: ' + op)

# --- Binning factor ---
root = tk.Tk()
root.withdraw()
proj_bin_size = simpledialog.askinteger(title="Projection binning", prompt="How many projections do you want to average?")
print('Every ' + str(proj_bin_size) + ' projections will be averaged.')

# --- Scan type ---
root = tk.Tk()
root.withdraw()
scan_type = simpledialog.askstring(title="Scan type", prompt="'Full' or 'Half' acquisition?")
print('Scan type: ' + scan_type + ' acquisition.')

# --- Sesame-only inputs (pixel size and energy not in metadata) ---
if source == 'SESAME':
    root = tk.Tk()
    root.withdraw()
    pixel_size = simpledialog.askfloat(title="Pixel size", prompt="Specify pixel size in units of um")
    print('Pixel size: ' + str(pixel_size) + ' um')

    root = tk.Tk()
    root.withdraw()
    energy = float(simpledialog.askinteger(title="Energy", prompt="Specify beam energy in units of keV"))
    print('Beam energy: ' + str(energy) + ' keV')


###############################################################################
# computation #################################################################
###############################################################################

start_time = datetime.now()
print('FYI: Conversion may take a while. How long depends on the size of your data and available resources.')
print('Start = ' + start_time.strftime('%Hh:%Mm:%Ss, %m.%d.%Y'))


###############################################################################
# source-specific: extract flats, darks, projections, metadata ################
###############################################################################

if source == 'BAM':

    data = fbam['entry/data/data']
    n = np.shape(data)[0]
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Input data reference acquired')
    print('Length of full data: ' + str(n))

    # Timestamp metadata
    ndarraytimestamp = np.array(fbam['entry/instrument/NDAttributes/NDArrayTimeStamp'])
    dd_ndarraytimestep = np.array(np.diff(ndarraytimestamp))
    mode = stat.mode(dd_ndarraytimestep)

    # Flats
    print("Extracting flats...", end=" ", flush=True)
    flats1_stack = np.array(data[:20])
    flats2_stack = np.array(data[n - 19:])
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Flats extracted')

    # Darks (BAM doesn't take darks, create synthetic ones)
    darks_stack = np.zeros_like(flats1_stack) + 100
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Darks created')

    # Identify and remove accelerating projections
    if 'SAMPLE_W' in fbam['/entry/instrument/NDAttributes']:
        samplew = np.array(fbam['/entry/instrument/NDAttributes/SAMPLE_W'])
    elif 'CT_MICOS_W' in fbam['/entry/instrument/NDAttributes']:
        samplew = np.array(fbam['/entry/instrument/NDAttributes/CT_MICOS_W'])
    samplew_radio = samplew[len(flats1_stack):-len(flats2_stack)]
    scanspan = int(10 * (np.floor(np.max(samplew_radio) / 10)))
    accspan = abs((np.max(samplew) - np.min(samplew)) - scanspan)
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Acceleration radiographs identified')

    n_allradio = n - 40
    n_acc = n_allradio % 100
    n_truradio = n_allradio - n_acc
    print('Number of all projections (excluding flats): ' + str(n_allradio))
    print('Number of final projections (after removing acceleration): ' + str(n_truradio))
    shift2tru = 20 + n_acc

    # Pixel size, energy, detector distance from metadata
    pixel_size_m  = fbam['entry/instrument/NDAttributes/CT_Pixelsize'][()][0] * 1e-6  # um -> m
    energy_val    = fbam['entry/instrument/NDAttributes/DMM_Energy'][()][0]           # keV
    det_dist_m    = (fbam['entry/instrument/NDAttributes/CT-Kamera-Z'][()][0] + 25) * 1e-4  # mm -> m

else:  # SESAME

    # Keep as lazy h5py dataset — same as BAM, no full load into RAM
    data = fbam['exchange/data']
    # First frame is a duplicate in Sesame files, skip it via shift
    shift2tru = 1
    n_truradio = data.shape[0] - 1  # exclude the duplicate first frame

    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Input data reference acquired')
    print('Number of projections: ' + str(n_truradio))

    # Flats and darks are stored separately in Sesame files (small, safe to load)
    print("Extracting flats and darks...", end=" ", flush=True)
    flats1_stack = np.array(fbam['exchange/data_white'])
    flats2_stack = np.array(fbam['exchange/data_white'])  # only one set; reuse for second slot
    darks_stack  = np.array(fbam['exchange/data_dark'])
    ctime = datetime.now()
    print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Flats and darks extracted')

    # Rotation angles from metadata
    angles_raw = np.array(fbam['exchange/theta'])[1:]
    scanspan = int(np.max(angles_raw))

    # Pixel size, energy, detector distance from user input
    pixel_size_m = pixel_size * 1e-6                                                   # um -> m
    energy_val   = energy                                                               # keV
    det_dist_m   = fbam['measurement/instrument/detector_motor_stack/detector_z'][0] * 1e-4  # mm -> m


###############################################################################
# shared: binning, section sizes, metadata arrays #############################
###############################################################################

if n_truradio % proj_bin_size == 0:
    no_bins = n_truradio // proj_bin_size
else:
    no_bins = int(n_truradio / proj_bin_size) + 1

half = no_bins // 2

# Alignment frame: midpoint projection (one chunk, read now)
# shift2tru accounts for skipped frames (acceleration in BAM, duplicate first frame in Sesame)
align_src_start = shift2tru + half * proj_bin_size
align_src_end   = min(align_src_start + proj_bin_size, shift2tru + n_truradio)
alignment_frame = np.mean(np.array(data[align_src_start:align_src_end]), axis=0)
alignment_stack = alignment_frame.reshape(1, x, y)

# Section sizes
n_darks     = len(darks_stack)
n_flats1    = len(flats1_stack)
n_radiosA   = half
n_radiosB   = no_bins - half
n_flats2    = len(flats2_stack)
n_alignment = 1
total_frames = n_darks + n_flats1 + n_radiosA + n_flats1 + n_radiosB + n_flats2 + n_alignment
print('Total frames in output: ' + str(total_frames))
print('First half of radios: '   + str(n_radiosA))
print('Second half of radios: '  + str(n_radiosB))

# Image key control
image_key_control = np.concatenate([
    [ImageKey.DARK_FIELD]  * n_darks,
    [ImageKey.FLAT_FIELD]  * n_flats1,
    [ImageKey.PROJECTION]  * n_radiosA,
    [ImageKey.FLAT_FIELD]  * n_flats1,
    [ImageKey.PROJECTION]  * n_radiosB,
    [ImageKey.FLAT_FIELD]  * n_flats2,
    [ImageKey.ALIGNMENT]   * n_alignment,
])
assert len(image_key_control) == total_frames
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Image key control created')
print('len image key control: ' + str(len(image_key_control)))

# Rotation angles
angles = np.linspace(0, scanspan, no_bins)
rotation_angle = np.concatenate([
    [0.0]          * n_darks,
    [0.0]          * n_flats1,
    angles[:n_radiosA],
    [scanspan / 2] * n_flats1,
    angles[n_radiosA:],
    [scanspan]     * n_flats2,
    [scanspan / 2] * n_alignment,
])
assert len(rotation_angle) == total_frames
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Rotation angles specified')
print('len of rotation angles: ' + str(len(rotation_angle)))


###############################################################################
# build NXtomo and save metadata ##############################################
###############################################################################

ureg = pint.UnitRegistry()

new_nx = NXtomo()
new_nx.instrument.detector.image_key_control = image_key_control
new_nx.sample.rotation_angle = rotation_angle * ureg.degree
new_nx.instrument.detector.field_of_view = str(scan_type)

# nxtomo >= 3.0 accepts pint quantities; older versions expect plain floats
_nxtomo_version = tuple(int(x) for x in pkg_version('nxtomo').split('.')[:2])
print('nxtomo version: ' + pkg_version('nxtomo'))
if _nxtomo_version >= (3, 0):
    new_nx.instrument.detector.x_pixel_size = pixel_size_m * ureg.meter
    new_nx.instrument.detector.y_pixel_size = pixel_size_m * ureg.meter
    new_nx.energy = energy_val * ureg.keV
    new_nx.instrument.detector.distance = det_dist_m * ureg.meter
else:
    new_nx.instrument.detector.x_pixel_size = float(pixel_size_m)
    new_nx.instrument.detector.y_pixel_size = float(pixel_size_m)
    new_nx.energy = float(energy_val)
    new_nx.instrument.detector.distance = float(det_dist_m)

ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Metadata specified')

# Output path
root_out_path = op + '/' + filename + '/'
if not os.path.exists(root_out_path):
    os.mkdir(root_out_path)
out_path = root_out_path + filename[:-3] + '_bin' + str(proj_bin_size) + '_v2_nxtomo.nx'
if os.path.exists(out_path):
    out_path = out_path[:-3] + '_1.nx'

new_nx.save(file_path=out_path, data_path='entry', overwrite=True)
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': NXtomo metadata saved. Now creating data dataset and writing chunk by chunk...')


###############################################################################
# create data dataset and write chunk by chunk ################################
###############################################################################

input_dtype = data.dtype

with h5py.File(out_path, 'r+') as f_out:

    ds = f_out.create_dataset(
        'entry/instrument/detector/data',
        shape=(total_frames, x, y),
        dtype=input_dtype,
        chunks=(1, x, y),
    )

    # Hard link expected by NXtomo/nabu
    f_out['entry/data/data'] = ds

    write_idx = 0

    # Darks
    for frame in tqdm(darks_stack, desc='Writing darks'):
        ds[write_idx] = frame
        write_idx += 1

    # Flats1
    for frame in tqdm(flats1_stack, desc='Writing flats1'):
        ds[write_idx] = frame
        write_idx += 1

    # RadiosA (first half of binned projections)
    for i in tqdm(range(n_radiosA), desc='Writing radiosA'):
        a = shift2tru + i * proj_bin_size
        b = min(a + proj_bin_size, shift2tru + n_truradio)
        chunk = np.array(data[a:b], dtype=np.float32)
        ds[write_idx] = np.mean(chunk, axis=0) if proj_bin_size > 1 else chunk[0]
        write_idx += 1

    # Flats1 again (between halves)
    for frame in tqdm(flats1_stack, desc='Writing flats1 (mid)'):
        ds[write_idx] = frame
        write_idx += 1

    # RadiosB (second half of binned projections)
    for i in tqdm(range(n_radiosA, no_bins), desc='Writing radiosB'):
        a = shift2tru + i * proj_bin_size
        b = min(a + proj_bin_size, shift2tru + n_truradio)
        chunk = np.array(data[a:b], dtype=np.float32)
        ds[write_idx] = np.mean(chunk, axis=0) if proj_bin_size > 1 else chunk[0]
        write_idx += 1

    # Flats2
    for frame in tqdm(flats2_stack, desc='Writing flats2'):
        ds[write_idx] = frame
        write_idx += 1

    # Alignment
    ds[write_idx] = alignment_stack[0]
    write_idx += 1

    assert write_idx == total_frames, f"Frame count mismatch: wrote {write_idx}, expected {total_frames}"

fbam.close()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': All data written into output file')

end_time = datetime.now()
print('End = ' + end_time.strftime('%Hh:%Mm:%Ss, %m.%d.%Y'))
duration = end_time - start_time
total_seconds = int(duration.total_seconds())
hours, remainder = divmod(total_seconds, 3600)
minutes, seconds = divmod(remainder, 60)
print(f'Nx conversion COMPLETE. Conversion duration = {hours}h:{minutes}m:{seconds}s')
print('Find your converted file in ' + out_path)
input('\nPress Enter to close...')
