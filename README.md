# MeerKAT-Cookbook
Recipes for MeerKAT data interaction and processing
* Accessing the data via the MeerKAT data archive
* Inspecting the data using the MeerKAT `katdal` software interface
* CASA recipes for processing MeerKAT data


## ARCHIVE
All MeerKAT data is accessed via the [SARAO archive](https://archive.sarao.ac.za/)     
MeerKAT data is stored in a flexible format called MeerKAT Visibility Format (MVF), and data is accessed and processed as needed.    
User guideline to register, access and retrieve data from the archive are provided through the      
[Archive Interface User Guide](https://archive.sarao.ac.za/statics/Archive_Interface_User_Guide.pdf)

**Note to reader**    
MeerKAT data files are large and combining the data for an observation using the full array into single files cause sizes of Giga- to Tera bytes.
These files are to big for standard io-operations.     

Some introduction usage example notebooks can be found in 
[archive](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/archive)

**TODO**
* Add workbook for MS creation from archive as soon as option is completed and tested.
Contrast with commandline instruction (hopefully with calibration option) shown in CASA recipe section


## KATDAL
Interacting with any of the MeerKAT observation files is made easy by the `katdal` python library   
[katdal repository](https://github.com/ska-sa/katdal)
Open source library available from PyPI
```
pip search katdal
katdal (0.13)  - Karoo Array Telescope data access library for interacting with data sets in the MeerKAT Visibility Format (MVF)
```
Detail `katdal` documentation with user guide instructions can be found on the [katdal's documentation](https://katdal.readthedocs.io/en/latest/index.html) page

`katdal` is specifically developed to allow efficient access to MeerKAT Visibility Format (MVF). It is fully integrated to access data via the `kat archive` line, optimised for large file data access and memory usage. 
Introductory Jupyter notebooks illustrating some example data interaction and inspection using 
[katdal](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/katdal)


## CASA
`katdal` provides a script to convert these data sets to CASA MeasurementSets.


 -fin-
