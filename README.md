# MeerKAT-Cookbook
Recipes for MeerKAT data interaction and processing


## ARCHIVE
MeerKAT data is stored in a flexible format called MeerKAT Visibility Format (MVF), and accessed/processed as needed.    
All MeerKAT data is accessed via the [SARAO archive](https://archive.sarao.ac.za/)     

**Important Notes**:
* MeerKAT archive is access restricted, requiring registration and login
* To access data in the MeerKAT archive a **token** is required.

User guideline to register, access and retrieve data from the archive are provided in the
[Archive Interface User Guide](https://archive.sarao.ac.za/statics/Archive_Interface_User_Guide.pdf)

Some introduction usage example notebooks can be found in 
[archive](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/archive)

**Note to reader**    
MeerKAT data files are large and combining the data for an observation using the full array into single files cause sizes of Giga- to Tera bytes.
These files are to big for standard io-operations.     

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


## UTILS
An assortment of helper scripts to assist the astronomer

 -fin-
