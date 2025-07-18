# Tomography data conversion tools

These scripts convert data collected at various synchrotron beamlines to data in .nx format, which is compatible with ESRF's tomography reconstruction software Tomwer (https://tomotools.gitlab-pages.esrf.fr/tomwer). \
These scripts create .nx files following guidelines described here: https://tomotools.gitlab-pages.esrf.fr/nxtomo/tutorials/create_from_scratch.html \
Currently there are two scripts, supporting conversion of data either from Bessy's BAMline or Sesame's BEATS beamline. 

### How to use (linux) ###
1. Download the .py file that you need (depending on where your data was collected, Bessy or Sesame) and the requirements.txt file, to a path of your choice.
2. Navigate to this path.
3. Install required packages: pip install -r requirements.txt 
4. Run the script: python ******.py 


Developed in the Zaslansky lab of the Department of Operative, Preventive and Pediatric Dentistry at Charité – Universitätsmedizin Berlin, with funding from the DFG (FOR5657). 
