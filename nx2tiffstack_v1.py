#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 09:44:07 2025

@author: adminisa
"""

# VERSION 1

# Developed in 08.2025 by Isabela Vitienes, in the Zaslansky lab 
# of the Department of Operative, Preventive and Pediatric Dentistry 
# at Charité – Universitätsmedizin Berlin, with funding from the 
# DFG (FOR5657).


##################################################################
# This code converts holotomo aligned projections in .nx format  #
# (pynx output) into a tiff-stack that can be reconstructed with # 
# alternative methods.                                           #
##################################################################

#############################################################
#                                                           #  
# INPUT: (1) User-selected .nx file of aligned projections  #
#        (2) User-selected output directory                 #
# OUTPUT: Folder in output directory containing tiff stack. #
#                                                           #
#############################################################



import numpy as np
import os
import tkinter as tk
from tkinter import filedialog
import time
import sys
import itertools
import tifffile as tif
from datetime import datetime
import h5py




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
input_path = filedialog.askopenfilename(title="Select an .nx file", filetypes=[("nx file", "*.nx")])
filename = upstream_path(input_path)[1][:-3]
fnx = h5py.File(input_path, 'r')


# Prompt user to select output directory
root = tk.Tk()
root.withdraw()  # Hide the root window
op = filedialog.askdirectory(title="Select output location")
op_f = op + '/' + filename
if not os.path.isdir(op_f):
    os.mkdir(op_f)
print('The output tiff stack will be stored here: ' + op_f)



start_time = datetime.now()
print('Loading projectins from input')
print('Start = ' + start_time.strftime('%Hh:%Mm:%Ss, %m.%d.%Y'))
data = np.array(fnx['/entry_1/data/data'])
z, x, y = data.shape
print('Projections loaded.')
print('Your dataset has ' + str(z-1) + ' projections.')
print('Each projection is ' + str(x) + ' x ' + str(y) + ' pixels.')


start_time = datetime.now()
print('Saving tiff slices')
print('Start = ' + start_time.strftime('%Hh:%Mm:%Ss, %m.%d.%Y'))
for i in range(z):
    
    if i<10:
        slcno = '000' + str(i)
    elif i<100:
        slcno = '00' + str(i)
    elif i<1000:
        slcno = '0' + str(i)
    else:
        slcno = str(i)
    
    path = op_f + '/' + filename + '_' + slcno + '.tif'
    
    tif.imwrite(path, data[i].astype(np.uint8))


print('Process complete.')
end_time = datetime.now()
print('End = ' + end_time.strftime('%Hh:%Mm:%Ss, %m.%d.%Y'))

