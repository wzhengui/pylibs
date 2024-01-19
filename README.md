# pylibs
* **Repository of python functions/classes/scripts. This library is designed to handle most of our routine work. Processing SCHISM related work is a significant part of the usgae of this library.**   

* In pylibs, there are mainly two groups of functions: <br>
  *  generic function: lpfilt, ReadNC, inside_polygon, proj,  ... <br>
  *  schism-related function:  sms2grd, read_schism_hgrid, read_schism_output), ... <br>

* Installation <br>
  * `pip install pylibs-ocean` (user mode)
  * `pip install pylibs-ocean[mpi,shapefile,projection,eof]` (comprehensive user mode)
  * `git clone https://github.com/wzhengui/pylibs.git; cd pylibs; pip install -e .` (developer mode)

* Usage <br>
  * explicit import:  `from pylib import zdata, ReadNC, read_schism_hgrid, sms2grd`
  * import mport:   `from pylib import *` (import all)

* Directories  <br>
  * Scripts: sample scripts for using pylibs <br>
  * Utility: python library functions <br>
    * pylib.py: tool for importing all necessary and frequently-used python functions/packages <br>
    * mylib.py: defined functions/classes  
    * schism_file.py: schism related functions/classes
