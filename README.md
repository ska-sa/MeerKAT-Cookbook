# MeerKAT-Cookbook
Recipes for MeerKAT data interaction and processing presented in CASA-Jupyter notebook format for easy
interaction.    
Usage and installation instruction for CASA-Jupyter can be found in [Github](https://github.com/aardk/jupyter-casa)

MeerKAT data is stored in a flexible format called MeerKAT Visibility Format (MVF), and accessed/processed as needed.    
Easy access and software specifically developed to handle MeerKAT large data sizes are provided through the MeerKAT archive and `katdal` packages.

**Important Notes**:    
If you are using the docker container to run the CASA-Jupyter notebooks, you will need to install
`katdal`
* run/start the docker container
* enter the docker container as root: docker exec -tiu root `container ID` bash
* install `katdal`: pip install katdal
* for active notebooks, restart the notebook kernel (using the `Kernel` menu option)
* running a cell with `import katdal` should now work

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
[Convert MVF dataset(s) to MeasurementSet](https://github.com/ska-sa/MeerKAT-Cookbook/blob/master/archive/Convert%20MVF%20dataset(s)%20to%20MeasurementSet.ipynb)

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

Positional astronomy calculation use the [PyEphem](https://pypi.org/project/ephem/) library

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


 -fin-
