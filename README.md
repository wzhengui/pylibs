# pylibs
* Repository of python functions. There are two types of functions: <br>
  *  For general purpose use (e.g., lpfilt, inside_polygon,proj,get_subplot_position) <br>
  *  For procesing SCHISM-Grid related input/outputs/analysis (e.g., read_schism_hgrid,read_schism_bpfile) <br>

* Directories  <br>
  * Scripts: sample scripts for using pylibs <br>
  * Utility: python library functions <br>
    * pylib.py: tool for importing all necessary python funcitons for daily routine work (not recommend for jobs with high efficiency demand) <br>
    * mylib.py: self-defined functions/classes  
    * schism_file.py: schism-grid related functions/classes

* Usage: <br>
  * 1). First, add pylibs path to environmental variable 'PYTHONPATH' <br>
      * e.g. on c-shell: setenv PYTHONPATH 'mydir/pylibs/Utility/:mydir/pylibs/Scripts/' <br>
  * 2). Second add following two lines in the beginning of each script  <br>
      * #!/usr/bin/evn python3 <br> from pylib import *  <br>


