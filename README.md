# MeerKAT-Cookbook
Recipes for MeerKAT data interaction and processing

MeerKAT data is stored in a flexible format called MeerKAT Visibility Format (MVF), and accessed/processed as needed.    
Easy access and software specifically developed to handle MeerKAT large data sizes are provided through the MeerKAT archive and `katdal` packages.

**Note to reader**    
MeerKAT data files are large and combining the data for an observation using the full array into single files cause sizes of Giga- to Tera bytes.
These files are to big for standard io-operations.     


## ARCHIVE
All MeerKAT data is accessed via the [SARAO archive](https://archive.sarao.ac.za/)     
**Important Notes**:
* MeerKAT archive is access restricted, requiring registration and login
* To access data in the MeerKAT archive a **token** is required.

User guideline to register, access and retrieve data from the archive are provided in the
[Archive Interface User Guide](https://archive.sarao.ac.za/statics/Archive_Interface_User_Guide.pdf)

Example notebooks showing data interaction and extraction methods can be found in the 
[archive](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/archive) folder
* Using tokens for `katdal` processing:
[Accessing MeerKAT observation data](https://github.com/ska-sa/MeerKAT-Cookbook/blob/master/archive/Accessing%20MeerKAT%20observation%20data.ipynb)
* `katdal` provides a script to convert these data sets to CASA MeasurementSets. 
Using tokens to convert MVF files to CASA MeasurementSet:
[Convert MVF dataset(s) to MeasurementSet](https://github.com/ska-sa/MeerKAT-Cookbook/blob/master/utils/Convert%20MVF%20dataset(s)%20to%20MeasurementSet.ipynb)

If you are following standard interferometric imaging data reduction using CASA measurement sets, you can also use the Direct Download Link to create and download a measurement set instead.    
See [Archive Interface User Guide](https://archive.sarao.ac.za/statics/Archive_Interface_User_Guide.pdf) for detail.


## KATDAL
Interacting with any of the MeerKAT observation files is made easy by the `katdal` python library   
[katdal repository](https://github.com/ska-sa/katdal)
Open source library available from PyPI
```
pip search katdal
katdal (0.13)  - Karoo Array Telescope data access library for interacting with data sets in the MeerKAT Visibility Format (MVF)
```    
Easy install     
```
pip install katdal
```

Detail `katdal` documentation with user guide instructions can be found on the `katdal`
[User guide](https://katdal.readthedocs.io/en/latest/index.html)

`katdal` is specifically developed to allow efficient access to MeerKAT Visibility Format (MVF).
It is fully integrated to access data via the `katarchive` line, optimised for large file data access and memory usage. 
Introductory Jupyter notebooks illustrating some example data interaction and inspection using 
[katdal](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/katdal)

Cutting edge functionality can be obtained by installing `katdal` directly from the GitHub repository    
```
pip install git+https://github.com/ska-sa/katdal.git
```
**Import Note**
Care should be taken since installing the master from GitHub might not be as stable as PyPI.


## CASA
CASA MeasurementSet data tables can be created using a convenient helper script `mvftoms.py` available from `katdal` installation.     
Measurement sets can be downloaded directly from the MeerKAT archive using some sensible defaults when created.     
Examples on how to create measurement sets from a user control environment using tokens from the archive are given in example notebooks in the utils folder.

Standard recipes for flagging, calibration and imaging are provided in the 
[casa](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/casa) folder.

The reader is advised to also consult the MeerKATHI pipeline.
When comparing the recipe to the MeerKATHI configuration file, close .... of the recipe will be observed.


## UTILS
An assortment of helper scripts to assist the astronomer

 -fin-
