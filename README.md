# Tomography data conversion tools

The  *h5tonx* script converts data collected at various synchrotron beamlines to data in .nx format, which is compatible with ESRF's tomography reconstruction software Tomwer (https://tomotools.gitlab-pages.esrf.fr/tomwer). \
The created .nx output file follows structure guidelines described here: https://tomotools.gitlab-pages.esrf.fr/nxtomo/tutorials/create_from_scratch.html \
Currently the script supports conversion of data either from Bessy's BAMline (pre and post 2022) or Sesame's BEATS beamline.

### System requirements ###
It can be run on minimal Linux or Windows systems.
It requires Python, and has been tested using Python 3.12 and 3.13.
File conversion requires $$R_{total} = \frac{(60 + n_{bin}) \cdot x \cdot y \cdot \frac{b}{8}}{1024^2} + 500 \text{ MB}$$ MBs of RAM, where \
- $x, y$ = slice dimensions in pixels
- $b$ = bit depth (8, 16, or 32)
- $n_{bin}$ = binning factor (number of projections averaged per chunk) \
It is strongly recommended to set up and execute within a conda environment.

### How to use (linux) ###
1. Download the *h5tonx.py* and *requirements.txt* file, to a path of your choice and navigate to this path.
3. Install required packages: *pip install -r requirements.txt* or *conda install -r requirements.txt*
4. Execute the python script: *python h5tonx.py* 

### How to use (windows) ###
1. Download *h5tonx.py*, *requirements.txt*, and *run_h5tonx.bat* to the same path, of your choice.
2. Run (double click) the *run_h5tonx.bat* file. This will automatically check if python is installed, if needed install packages listed in *requirements.txt*, and run *h5tonx.py*.\
*Note: It is possible to set a proxy in the .bat file if needed.*

---

First developed in 06.08.2025 by Isabela Vitienes, in the Zaslansky lab of the Department of Operative, Preventive and Pediatric Dentistry at Charité – Universitätsmedizin Berlin, with funding from the DFG (FOR5657). 
Last updated: 13.03.2026
