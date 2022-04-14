# pylibs
* **Repository of python functions/classes/scripts. This library is designed for dealing with different kinds of daily work. Processing SCHISM related work is a significant part of the usgae of this library.**   

* There are mainly two types of functions: <br>
  *  For general purpose use (e.g., lpfilt, inside_polygon,proj,get_subplot_position) <br>
  *  For pre/post-procesing SCHISM-Grid related input/outputs/analysis (e.g., read_schism_hgrid,read_schism_bpfile) <br>

* Directories  <br>
  * Scripts: sample scripts for using pylibs <br>
  * Utility: python library functions <br>
    * pylib.py: tool for importing all necessary python funcitons for daily routine work <br>
    * mylib.py: self-defined functions/classes  
    * schism_file.py: schism-grid related functions/classes

* Usage: <br>
  * Install pylibs: `cd mydir; git clone https://github.com/wzhengui/pylibs.git; pip install -e mydir/pylibs`
  * Import pylibs: 1). `from pylib import ReadNC, read_schism_hgrid ` (explicit), 2). `from pylib import *` (import all)


